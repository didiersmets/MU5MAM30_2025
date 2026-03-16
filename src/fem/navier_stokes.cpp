#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"
#include "boundary.h"
#include "P1.h"
#include "tiny_blas.h"
#include "mesh.h"

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
	, r(N)
	, p(N)
	, Ap(N)
	,is_bnd(N)
	,b(N)// right-hand side values at vertices

{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);
#endif
	compute_boundary_vertices(m,is_bnd);	

	inited = false;
	t = 0;
}


void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));// zero transport, T has size N, one value per vertex

	for (size_t t = 0; t < m.triangle_count(); t++) {// for each triangle
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		assert(a < N && b < N && c < N);//verify that the indices are valid
		double sum = omega[a] + omega[b] + omega[c];// sum of the vorticity at the vertices of the triangle
		T[a] += sum * (psi[b] - psi[c]);
		T[b] += sum * (psi[c] - psi[a]);
		T[c] += sum * (psi[a] - psi[b]);
	}
	for (size_t v = 0; v < N; v++) {
		T[v] *= 1.0 / 6;
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	double b2, r2, rel_error;
	size_t iter;

	double *R = r.data;
	double *Pv = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;
	double *Psi = psi.data;
	double *B = b.data;// right-hand side values at vertices


	
	// 0) we build S, we will modify it imposing boundary condition 
	#if !USE_FEM_MATRIX
	build_P1_stiffness_matrix(m, P,S);
	#endif
	//I rebuild S every time , otherwise after first temporal step the system is no longer the original one.
	
	M.mvp(Om, MOm);
	// 1) B = M * omega
	blas_copy(MOm, B, N);// B = M * omega;
	// access to matrix entries
	uint32_t *rs = S.row_start;
	uint32_t *cc = S.col;
	double   *Sd = S.data.data;
	
	// 2) zero columns corresponding to boundary vertices, on internal rows
	for (uint32_t i = 0; i < (uint32_t)N; ++i) {
		if (is_bnd[i]) continue;

		for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
			uint32_t j = cc[k];
			if (is_bnd[j])
				Sd[k] = 0.0;
		}
	}
	// 3) set boundary rows equal to identity
	for (uint32_t i = 0; i < (uint32_t)N; ++i) {
		if (!is_bnd[i]) continue;


		for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
			if (cc[k] == i) {
				Sd[k] = 1.0;
			}
			else {
				Sd[k] = 0.0;
			}
		}

		B[i] = 0.0;
	}


	/* Compute rhs norm2 */
	b2 = blas_dot(B, B, N);

	/* Form initial R and P */
	S.mvp(Psi, R);
	blas_axpby(1, B, -1, R, N);
	blas_copy(R, Pv, N);
	r2 = blas_dot(R, R, N);
	rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);

	/* Iterate until convergence */
	iter = 0;
	while ((rel_error > tol) && (iter++ < iter_max)) {

		/* Compute AP */
		S.mvp(Pv, AP);

		/* Update Psi */
		double alpha = r2 / blas_dot(Pv, AP, N);
		blas_axpy(alpha, Pv, Psi, N);
		/* Update R */
		blas_axpy(-alpha, AP, R, N);

		/* Update r2 and P */
		double beta = 1.0 / r2;
		r2 = blas_dot(R, R, N);
		rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);
		beta *= r2;
		blas_axpby(1, R, beta, Pv, N);
	}

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	double b2, r2, rel_error;

	size_t iter1, iter2;

	double *R = r.data;
	double *Pv = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;

	iter1 = compute_stream_function();
	#if !USE_FEM_MATRIX
		build_P1_stiffness_matrix(m, P, S);
	#endif
	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Form rhs, saved in P */
	compute_transport(Pv);
	
	blas_axpby(1, MOm, dt, Pv, N); //Pi = 1 * Momi + dt * Pi M*omega(t) + dt * T(Omega,Psi)(t) on first N values of P
	b2 = blas_dot(Pv, Pv, N);//compute the norm of the rhs, used for convergence check, on first N values of P

	/* Form initial R and P */
	S.mvp(Om, R);// R = S * omega(t)
	blas_axpby(1, MOm, dt * nu, R, N);// R =  M * omega(t)+ dt * nu * S * omega(t) = (M + \nu * dt * S)omega(t)
	blas_axpby(1, Pv, -1, R, N);// R = Pv - R = b-A*omega(t) since in Pv I've the exct rhs, R is the initial residual
	blas_copy(R, Pv, N);// P = R, initial search direction is the residual
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	/* Iterate until convergence (and at least once) */
	iter2 = 0;
	do {

		/* Compute AP (invalidates Mom) */
		S.mvp(Pv, AP);
		M.mvp(Pv, MOm); /* MOm used as temp storage */
		blas_axpby(1, MOm, dt * nu, AP, N);

		/* Update Om */
		double alpha = r2 / blas_dot(Pv, AP, N);
		blas_axpy(alpha, Pv, Om, N);

		/* Update R */
		blas_axpy(-alpha, AP, R, N);

		/* Update r2 and P */
		double beta = 1.0 / r2;
		r2 = blas_dot(R, R, N);
		rel_error = sqrt(r2 / b2);
		beta *= r2;
		blas_axpby(1, R, beta, Pv, N);

		/* Update MOm */
		M.mvp(Om, MOm);

		iter2++;
	} while ((rel_error > tol) && (iter2 <= iter_max));

	

	t += dt;

	(void)iter1;
	//printf("Iter 1 : %zu, Iter2 : %zu\n", iter1, iter2);
}
