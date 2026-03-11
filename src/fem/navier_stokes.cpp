#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "conjugate_gradient.h"

#include "P1.h"
#include "P2.h"
#include "tiny_blas.h"

#include "quadrature.h"

NavierStokesSolver::NavierStokesSolver(Mesh &m, const bool &use_fem_P2)
	: m(m), omega(0), Momega(0), psi(0), r(0), p(0), Ap(0), use_fem_P2(use_fem_P2)
{
#if USE_FEM_MATRIX
	if (use_fem_P2)
		throw std::runtime_error("P2 Lagrange elements are not available for FEM matrix framework");
	else
	{
		build_P1_mass_matrix(m, M);
		build_P1_stiffness_matrix(m, S);
	}
#else
	if (use_fem_P2)
	{

		m.build_edges();
		build_P2_CSRPattern(m, P);
		build_P2_mass_matrix(m, P, M);
		build_P2_stiffness_matrix(m, P, S);

		precompute_phi_at_quad_points(phi_at_quad);
		precompute_grad_phi_at_quad_points(grad_phi_at_quad);
	}
	else
	{
		build_P1_CSRPattern(m, P);
		build_P1_mass_matrix(m, P, M);
		build_P1_stiffness_matrix(m, P, S);
	}
#endif
	rel_error = (double *)malloc(sizeof(double));
	t = 0;
	vol = M.sum();
	inited = false;

	N = M.cols;
	omega.resize(N);
	Momega.resize(N);
	psi.resize(N);
	r.resize(N);
	p.resize(N);
	Ap.resize(N);
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	/* We use the formula : \bar v = (\sum_{i, j} V_i * M_{ij})  / vol
	Then we set V <- V - \bar v */
	double *TEMP = (double *)malloc(N * sizeof(double));
	M.mvp(V, TEMP);
	double sum = blas_sum_in_place(TEMP, N);
	for (size_t i = 0; i < N; i++)
		V[i] -= sum / vol;
}

void NavierStokesSolver::compute_transport(double *T)
{
	throw std::runtime_error("This function should not be called.");
}

void NavierStokesSolver::compute_transport_P1(double *T)
{
	/* We use the formula :
	\forall j \in I, T[j] = \sum_{i, k} \Omega_i * \Psi_k \int_{\Omega} \phi_i * (\nabla^T \phi_k . \nabla \phi_j) */
	memset(T, 0, N * sizeof(double));

	double *OMEGA = omega.data;
	double *PSI = psi.data;

	size_t nt = m.triangle_count();
	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = m.indices[3 * tri];
		uint32_t b = m.indices[3 * tri + 1];
		uint32_t c = m.indices[3 * tri + 2];

		double omega_sum = OMEGA[a] + OMEGA[b] + OMEGA[c];

		/* Compute the explicit contribution for each node of the triangle */
		T[a] += (omega_sum * (PSI[c] - PSI[b])) / 6;
		T[b] += (omega_sum * (PSI[a] - PSI[c])) / 6;
		T[c] += (omega_sum * (PSI[b] - PSI[a])) / 6;
	}
}

void NavierStokesSolver::precompute_phi_at_quad_points(double (&phi_at_quad)[7][6])
{
	for (int q = 0; q < 7; ++q)
	{
		double L1 = Quadrature::lambda1[q];
		double L2 = Quadrature::lambda2[q];
		double L3 = Quadrature::lambda3[q];

		// Vertices
		phi_at_quad[q][0] = L1 * (2.0 * L1 - 1.0);
		phi_at_quad[q][1] = L2 * (2.0 * L2 - 1.0);
		phi_at_quad[q][2] = L3 * (2.0 * L3 - 1.0);

		// Edges
		phi_at_quad[q][3] = 4.0 * L1 * L2;
		phi_at_quad[q][4] = 4.0 * L2 * L3;
		phi_at_quad[q][5] = 4.0 * L3 * L1;
	}
}

void NavierStokesSolver::precompute_grad_phi_at_quad_points(Vec2d (&grad_phi_at_quad)[7][6])
{
	for (int q = 0; q < 7; ++q)
	{
		double L1 = Quadrature::lambda1[q];
		double L2 = Quadrature::lambda2[q];
		double L3 = Quadrature::lambda3[q];

		// Vertices
		grad_phi_at_quad[q][0] = {4.0 * L1 - 1, 0};
		grad_phi_at_quad[q][1] = {0, 4.0 * L2 - 1};
		grad_phi_at_quad[q][2] = {-(4.0 * L3 - 1.0), -(4.0 * L3 - 1.0)};

		// Edges
		grad_phi_at_quad[q][3] = {4.0 * L2, 4.0 * L1};
		grad_phi_at_quad[q][4] = {-4.0 * L2, 4 * (L3 - L2)};
		grad_phi_at_quad[q][5] = {4 * (L3 - L1), -4.0 * L1};
	}
}

void NavierStokesSolver::compute_transport_P2(double *T)
{
	/* We use the formula :
	\forall j \in I, T[j] = \sum_{i, k} \Omega_i * \Psi_k \int_{\Omega} \phi_i * (\nabla^T \phi_k . \nabla \phi_j) */

	memset(T, 0, N * sizeof(double));

	double *OMEGA = omega.data;
	double *PSI = psi.data;

	double w[7];
	for (size_t q = 0; q < 7; q++)
		w[q] = Quadrature::weights[q];

	size_t nt = m.triangle_count();
	size_t nv = m.vertex_count();
	for (size_t tri = 0; tri < nt; tri++)
	{
		/* Storage of the six DOFs of the current triangle */
		uint32_t nodes[6];
		nodes[0] = m.indices[3 * tri];
		nodes[1] = m.indices[3 * tri + 1];
		nodes[2] = m.indices[3 * tri + 2];
		nodes[3] = nv + *m.edge_idx.get(Edge(nodes[0], nodes[1]));
		nodes[4] = nv + *m.edge_idx.get(Edge(nodes[1], nodes[2]));
		nodes[5] = nv + *m.edge_idx.get(Edge(nodes[2], nodes[0]));

		Vec3 x1 = m.positions[nodes[0]];
		Vec3 x2 = m.positions[nodes[1]];
		Vec3 x3 = m.positions[nodes[2]];

		/* Computation of the geometry of the triangle */
		Vec3d e1 = {(double)x1[0] - (double)x3[0],
					(double)x1[1] - (double)x3[1],
					(double)x1[2] - (double)x3[2]};
		Vec3d e2 = {(double)x2[0] - (double)x3[0],
					(double)x2[1] - (double)x3[1],
					(double)x2[2] - (double)x3[2]};

		double G11 = norm2(e1);
		double G12 = dot(e1, e2);
		double G22 = norm2(e2);

		Vec3d n = cross(e1, e2);
		double detG = norm2(n);
		double sqrt_detG = norm(n);
		normalized(n);

		for (size_t q = 0; q < 7; q++)
		{
			double A_K_at_q = 0.0;
			Vec3d B_K_at_q = {0.0, 0.0, 0.0};
			Vec3d grad_phi_real_at_quad[6];
			for (size_t i = 0; i < 6; i++)
			{
				/* Computation of A_K(\xi_q) */
				A_K_at_q += OMEGA[nodes[i]] * phi_at_quad[q][i];

				/* Computation of B_K(\xi_q) */
				double ux = (1 / detG) * (G22 * grad_phi_at_quad[q][i].x - G12 * grad_phi_at_quad[q][i].y);
				double uy = (1 / detG) * (-G12 * grad_phi_at_quad[q][i].x + G11 * grad_phi_at_quad[q][i].y);

				grad_phi_real_at_quad[i] = ux * e1 + uy * e2;

				B_K_at_q += PSI[nodes[i]] * grad_phi_real_at_quad[i];
			}

			Vec3d n_cross_B_K = cross(n, B_K_at_q);

			/* Assemble all the elements in T */
			for (size_t i = 0; i < 6; i++)
			{
				T[nodes[i]] += w[q] * A_K_at_q * dot(n_cross_B_K, grad_phi_real_at_quad[i]) * sqrt_detG;
			}
		}
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/**********************************************************************
	 * Solve the system :
	 *
	 *  S * \Psi(t) = M * \Omega(t)
	 *
	 *********************************************************************/

	double *PSI = psi.data;
	double *MOMEGA = Momega.data;
	double *OMEGA = omega.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	M.mvp(OMEGA, MOMEGA);

	/* Solve using CG */
	iter = conjugate_gradient_solve(S, MOMEGA, PSI, R, P, AP, rel_error, tol, iter_max, false);

	set_zero_mean(PSI);

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	double *T = (double *)malloc(N * sizeof(double));
	compute_stream_function();
	if (use_fem_P2)
		compute_transport_P2(T);
	else
		compute_transport_P1(T);

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S) * omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	double *MOMEGA = Momega.data;
	double *OMEGA = omega.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	double *RHS = (double *)malloc(N * sizeof(double));
	
	blas_copy(MOMEGA, RHS, N);
	blas_axpy(dt, T, RHS, N);

	/* We directly solve the CG here to avoid creating the matrix M + \nu * dt * S */
	double *TEMP = (double *)malloc(N * sizeof(double));

	b2 = blas_dot(RHS, RHS, N);

	/* r_0 = RHS - M * omega(t) - \nu * dt * S * omega(t) */
	M.mvp(OMEGA, R);
	S.mvp(OMEGA, TEMP);
	blas_axpby(1, RHS, -1., R, N);
	blas_axpby(-nu * dt, TEMP, 1., R, N);

	/* p_0 = r_0 */
	blas_copy(R, P, N);

	r2 = blas_dot(R, R, N);
	*rel_error = sqrt(r2 / b2);

	size_t iter = 0;
	while ((iter < iter_max) && (*rel_error > tol))
	{
		/* Computation of (M + \nu * dt * S) * p_n */
		M.mvp(P, AP);
		S.mvp(P, TEMP);
		blas_axpby(nu * dt, TEMP, 1., AP, N);

		/* Computation of \alpha_n */
		double p_A2 = blas_dot(P, AP, N);
		double alpha = r2 / p_A2;

		/* Computation of omega_{n+1} */
		blas_axpby(alpha, P, 1., OMEGA, N);

		/* Computation of r_{n+1} */
		blas_axpby(-alpha, AP, 1., R, N);

		/* Computation of \beta_{n+1} */
		double new_r2 = blas_dot(R, R, N);
		double beta = new_r2 / r2;

		/* Computation of p_{n+1} */
		blas_axpby(1., R, beta, P, N);

		r2 = new_r2;
		*rel_error = sqrt(r2 / b2);
		iter++;
	}

	set_zero_mean(OMEGA);

	t += dt;
}
