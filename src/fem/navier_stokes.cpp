#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"
#include "conjugate_gradient.h"

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
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	/* We use the formula : \bar v = (\sum_{i, j} V_i * M_{ij})  / vol 
	Then we set V <- V - \bar v */
	double *TEMP = nullptr;
	M.mvp(V, TEMP);
	double sum = blas_sum_in_place(TEMP, N);
	for (size_t i = 0; i < N; i++)
		V[i] -= sum / vol;
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */
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

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	set_zero_mean(omega.data);

	t += dt;
}
