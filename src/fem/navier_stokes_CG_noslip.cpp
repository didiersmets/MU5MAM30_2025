#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <map>
#include <utility>
#include <cmath>
#include <vector>
#include <algorithm>

#include "navier_stokes_CG_noslip.h"
#include "P1.h"
#include "tiny_blas.h"
#include "projected_conjugate_gradient.h"
#include "matrix_sum.h"

//though i named it "non slip", this is non-penetrating + non-slip
//in $\Omega$ ：$\Delta\Psi = \omega$ (this always works, not special)
//   $\partial\Omega$ ：$\Psi = 0$ and $u = 0$ (they are special)
//so i want to change compute_stream_function and time_step, before this, i add how to distinguish border
//when there is no border like in sphere, use CG normal(set_zero_mean)

/*
* for omega, though we know u=0 on the border,omega=curl u is unknown, here i use Thom's Formula
* just like PCG for psi in penetrating, when we have omega_boundary_new, this becomes a non-homogene Diriclet problem,
* i use PCG again to ensure the omega_boundary_new remains the same while i calculate omega_inside_new
* luckly, i have already identified the boundary points that can be used directly
*/

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m), N(m.vertex_count()), omega(N), Momega(N), psi(N), r(N), p(N), Ap(N)
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);
#endif
	vol = M.sum();
	inited = false;
	t = 0;

	// identify the boundary
	is_boundary.assign(N, false);
	std::map<std::pair<uint32_t,uint32_t>, int> edge_count;
	for (size_t t = 0; t < m.triangle_count(); t++) {
		uint32_t a = m.indices[3*t+0], b = m.indices[3*t+1], c = m.indices[3*t+2];
		uint32_t tri[3][2] = {{a,b},{b,c},{c,a}};
		for (auto &e : tri) {
			uint32_t u = e[0] < e[1] ? e[0] : e[1];
			uint32_t v = e[0] < e[1] ? e[1] : e[0];
			edge_count[{u,v}]++;
		}
	}
	has_boundary = false;
	for (auto &kv : edge_count) {
		if (kv.second == 1) {
			is_boundary[kv.first.first]  = true;
			is_boundary[kv.first.second] = true;
			has_boundary = true;
		}
	}
}

void NavierStokesSolver::set_zero_mean(double *V)
{
    M.mvp(V, Ap.data);
    double s = 0.0;
    for (size_t i = 0; i < N; i++) s += Ap.data[i];
    double c = s / vol;
    for (size_t i = 0; i < N; i++) V[i] -= c;
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));
	for (size_t t = 0; t < m.triangle_count(); t++) {
		uint32_t a = m.indices[3 * t + 0], b = m.indices[3 * t + 1], c = m.indices[3 * t + 2];
		double sum = omega[a] + omega[b] + omega[c];
		T[a] += sum * (psi[c] - psi[b]);
		T[b] += sum * (psi[a] - psi[c]);
		T[c] += sum * (psi[b] - psi[a]);
	}
	for (size_t v = 0; v < N; v++) T[v] *= (1.0 / 6.0);
}

size_t NavierStokesSolver::compute_stream_function()
{
	M.mvp(omega.data, Momega.data);
	for (size_t i = 0; i < N; i++) Momega.data[i] = -Momega.data[i];

	if (!has_boundary) {
		set_zero_mean(Momega.data);
	} else {
        // same as in non-penetrating
		for (size_t i = 0; i < N; i++) {
			if (is_boundary[i]) {
				Momega.data[i] = 0.0;
				psi.data[i]    = 0.0;
			}
		}
	}

	double rel_error = 0.0;
	size_t iter = projected_conjugate_gradient_solve(
		S, Momega.data, psi.data, r.data, p.data, Ap.data,
		&rel_error, tol, iter_max, false,
        has_boundary, has_boundary ? &is_boundary : nullptr
	);

	if (!has_boundary) set_zero_mean(psi.data);
	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();

    // 2. non-slip: Thom's Formula
    if (has_boundary) {
        std::vector<double> omega_thom(N, 0.0);
        
        // Adjacency List
        std::vector<std::vector<uint32_t>> adj(N);
        for(size_t t = 0; t < m.triangle_count(); ++t) {
            uint32_t a = m.indices[3*t+0], b = m.indices[3*t+1], c = m.indices[3*t+2];
            adj[a].push_back(b); adj[a].push_back(c);
            adj[b].push_back(a); adj[b].push_back(c);
            adj[c].push_back(a); adj[c].push_back(b);
        }

        // use psi to calculate omega on the border
        for (size_t i = 0; i < N; ++i) {
            if (is_boundary[i]) {
                double sum_val = 0.0;
                int count = 0;
                
                // adjacency
                std::sort(adj[i].begin(), adj[i].end());
                adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());

                for (uint32_t nb : adj[i]) {
                    if (!is_boundary[nb]) { // focus on the inner points
                        double dx = m.positions[i][0] - m.positions[nb][0];
                        double dy = m.positions[i][1] - m.positions[nb][1];
                        double dz = m.positions[i][2] - m.positions[nb][2];
                        double dist2 = dx*dx + dy*dy + dz*dz;
                        
                        if (dist2 > 1e-12) {
                            sum_val += 2.0 * psi.data[nb] / dist2; // Thom's formula part1
                            count++;
                        }
                    }
                }
                if (count > 0) {
                    omega_thom[i] = sum_val / count; // Thom's formula part2
                }
            }
        }
        
        // Thom's formula part3
        for (size_t i = 0; i < N; ++i) {
            if (is_boundary[i]) omega.data[i] = omega_thom[i];
        }
    }

    // 3. combine RHS
    std::vector<double> T_vec(N, 0.0);
    std::vector<double> b_vec(N, 0.0);
    compute_transport(T_vec.data());                    
    M.mvp(omega.data, b_vec.data());                    
    blas_axpy(dt, T_vec.data(), b_vec.data(), N);       

    MatrixSum A(M, S, nu * dt, N); 
    double rel_error = 0.0;

    // 4. Non-homogeneous Dirichlet
    if (has_boundary) {
        std::vector<double> omega_bdry_only(N, 0.0);
        std::vector<double> A_omega_bdry(N, 0.0);
        
        for(size_t i = 0; i < N; ++i) {
            if(is_boundary[i]) omega_bdry_only[i] = omega.data[i];
        }
        
        //  A * omega_boundary
        A.mvp(omega_bdry_only.data(), A_omega_bdry.data());
        
        // b_I = b_I - A_IB * x_B
        for(size_t i = 0; i < N; ++i) {
            b_vec[i] -= A_omega_bdry[i];
        }

        // use Projection to limit the border
        for(size_t i = 0; i < N; ++i) {
            if(is_boundary[i]) b_vec[i] = 0.0;
        }
    }

    projected_conjugate_gradient_solve(
        A, b_vec.data(), omega.data, r.data, p.data, Ap.data,
        &rel_error, tol, iter_max, false,
        has_boundary, has_boundary ? &is_boundary : nullptr
    );

	if (!has_boundary) set_zero_mean(omega.data);
	t += dt;
}