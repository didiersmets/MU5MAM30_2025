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
	/* Your implementation goes here */

	/* Note that we only know how to compute \int{f_h * g_h} but here we can recover the mean 
	 * by computing \int{f_h * 1_h} where 1_h = \sum_{i\inV}{\phi_i} which is exactly the sum
	 * of the coefficients of M @ F where M is the mass matrix and F the vector induced by f_h
	 */

	/* To impose that V as mean zero we remove its mean from every component */
	M.mvp(V, Momega.data); // We use Momega as temp storage
	double sum = blas_sum_in_place(Momega.data, N); // /!\ destroy Momega.data
	double mean = sum / vol;
	
	for (size_t i=0; i<N; i++) {
		V[i] -= mean;
	}
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */

	/* We make use of the fact that T has a nice formula on the triangle ABC
	 * for x,y,z \in {A,B,C}, T =  (+/- with direct orientation) psi[x] * omega[y] * 1/6 * 1_{y!=z}
	 */

	/* For each triangle ijk in the mesh we add its contribution to 
	 * T[i] += (psi[i]*omega[j]-psi[j]*omega[i] + psi[j]*omega[k]-psi[k]*omega[j] + (!) psi[k]*omega[i] - (!) psi[i]*omega[k]) / 6;
	 * and subsequent contributions are given to other vertex of the triangle
	 */

	for (size_t tri_index=0; tri_index<m.index_count(); tri_index+=3) {
		uint32_t i = m.indices[tri_index];
		uint32_t j = m.indices[tri_index+1];
		uint32_t k = m.indices[tri_index+2];

		double tri_contribution = omega[i] + omega[j] + omega[k];

		T[i] += tri_contribution * (psi[j] - psi[k]);
		T[j] += tri_contribution * (psi[k] - psi[i]);
		T[k] += tri_contribution * (psi[i] - psi[j]);
	}

	for (size_t v = 0; v < N; v++) {
		T[v] *= 1.0/6;
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/* Your implementation goes here */

	/* We must solve S@psi = -M@omega 
	 * One easy solution is to use conjugate gradient here
	 * however pseudo-inverse is computationaly better as
	 * S is the same for every time-step
	 */

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

	/* Your implementation goes here */

	set_zero_mean(omega.data);

	t += dt;
}
