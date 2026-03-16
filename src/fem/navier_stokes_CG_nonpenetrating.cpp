#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <map>
#include <utility>
#include <cmath>
#include <vector>

#include "navier_stokes_CG_nonpenetrating.h"
#include "P1.h"
#include "tiny_blas.h"
#include "conjugate_gradient.h"   //for no boundary       
#include "projected_conjugate_gradient.h"  //when there is a boundary
#include "matrix_sum.h"                   // to combine M + nu*dt*S

//in $\Omega$ ：$\Delta\Psi = \omega$ (this always works, not special)
//   $\partial\Omega$ ：$\Psi = 0$ (this is special)
//so i want to change compute_stream_function on the border, before this, i add how to distinguish border
//when there is no border like in sphere, use CG normal(set_zero_mean)

/*
* how to add the BC is a problem, if i just let psi=0 after calculating, 
* it will be too sharp on the border, affecting points inside the mesh,
* so i use projected_conjugate_gradient_solve to avoid CG on boungary
*/

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
	, r(N)
	, p(N)
	, Ap(N)
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

	//identify the boundary points
	is_boundary.assign(N, false);
	std::map<std::pair<uint32_t,uint32_t>, int> edge_count;
	for (size_t t = 0; t < m.triangle_count(); t++) {
		uint32_t a = m.indices[3*t+0];
		uint32_t b = m.indices[3*t+1];
		uint32_t c = m.indices[3*t+2];
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

//for no boundary only
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
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
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
	for (size_t i = 0; i < N; i++)
		Momega.data[i] = -Momega.data[i];

	if (!has_boundary) {
		set_zero_mean(Momega.data);
	} else {
		/* to avoid machine error */
		for (size_t i = 0; i < N; i++) {
			if (is_boundary[i]) {
				Momega.data[i] = 0.0;
				psi.data[i]    = 0.0;
			}
		}
	}

	double rel_error = 0.0;
    // PCG here
	size_t iter = projected_conjugate_gradient_solve(
		S, Momega.data, psi.data, r.data, p.data, Ap.data,
		&rel_error, tol, iter_max, false,
        has_boundary, has_boundary ? &is_boundary : nullptr
	);

	//CG here
	if (!has_boundary) {
		set_zero_mean(psi.data);
	}

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();

    // combine RHS: b = M*omega + dt*T
    std::vector<double> T_vec(N, 0.0);
    std::vector<double> b_vec(N, 0.0);

    compute_transport(T_vec.data());                    
    M.mvp(omega.data, b_vec.data());                    
    blas_axpy(dt, T_vec.data(), b_vec.data(), N);       

    // combine A = M + nu*dt*S
    MatrixSum A(M, S, nu * dt, N); 
    double rel_error = 0.0;

    conjugate_gradient_solve(
        A, b_vec.data(), omega.data, r.data, p.data, Ap.data,
        &rel_error, tol, iter_max, false
    );

	if (!has_boundary) {
		set_zero_mean(omega.data);
    }

	t += dt;
}