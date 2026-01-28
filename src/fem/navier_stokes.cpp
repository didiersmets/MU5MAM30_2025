#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

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
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	M.mvp(V,Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i) {
		V[i] -= s / vol;
	}
}

void NavierStokesSolver::compute_transport(double *T,double dt)
{
	// computes T(omega,psi) and adds it to T
	// memset(T, 0, N * sizeof(double));
	for (size_t t=0;t<m.index_count();t+=3){
		uint32_t Ai = m.indices[t];
		uint32_t Bi = m.indices[t+1];
		uint32_t Ci = m.indices[t+2];
		uint32_t points[3] = {Ai,Bi,Ci};
		Vec3 AB = B-A;
		Vec3 AC = C-A;
		double sum = omega[Ai] + omega[Bi] + omega[Ci];
		for (uint32_t k =0;k<3;k++)
			T[points[k]]+= dt * sum * (psi[points[(k-1)%3]] - psi[points[(k+1)%3]]);
		}
}

size_t NavierStokesSolver::compute_stream_function()
{
	// solves for psi in : S psi = -M omega
	size_t iter = 0;
	M.mvp(omega,Momega)

	double rel_error;
	iter = conjugate_gradient_solve(
			S,
			Momega,
			psi,
			r,
			p,
			Ap
			&rel_error,
			10e-6,
			500,
			inited);
	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();
	compute_transport(Momega,dt)

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Your implementation goes here */

	set_zero_mean(omega.data);

	t += dt;
}
