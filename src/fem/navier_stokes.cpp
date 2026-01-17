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
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */
	double* v_omega = omega.data;
	double* v_psi = psi.data;

	/*

	T(Omega(t), Psi(t))_k = Sum_{i,j \in [1,N]} Omega_i(t) * Psi_j(t) * \integral psi_i * grad_perp(psi_j) \dot grad(psi_k) dA

	The two gradients are constants and the integral of psi is 1/6.
	T_A = 2 * A(T) / 6 * (\sum_{i in triangle} Omega_i) * (\sum_{j in triangle} Psi_j * grad_perp(psi_j) \dot grad(psi_A))
						 ^-- is computed once [sum_omega]						^						   ^-- is Zero for j = A, +1 or -1 in the other cases
																				^-- is stored in the v_psi array
	
																				The integral over the whole domain can be splitted as a sum of integrals over each triangle.
	We conpute this on each triangle and accumulate the results.

	The computation is done on a standard triangle, the transofrmation is up to a constant factor.

	T_i =
	i \in {A,B,C}
	
	
	*/

	// iterate over all triangles
	size_t tri_count = m.triangle_count();
	for (size_t t = 0; t < tri_count; t++) {
		// get vertex indices
		uint32_t ia = m.indices[3 * t + 0];
		uint32_t ib = m.indices[3 * t + 1];
		uint32_t ic = m.indices[3 * t + 2];

		// get vertex positions
		Vec3f A = m.positions[ia];
		Vec3f B = m.positions[ib];
		Vec3f C = m.positions[ic];

		// compute edges
		Vec3f AB = B - A;
		Vec3f AC = C - A;

		// area of the triangle
		double area = 0.5 * norm(cross(AB, AC));

		double sum_omega = v_omega[ia] + v_omega[ib] + v_omega[ic];
		// compute local contributions to transport term
		T[ia] += (area / 3.0) * sum_omega * (v_psi[ic] - v_psi[ib]);
		T[ib] += (area / 3.0) * sum_omega * (v_psi[ia] - v_psi[ic]);
		T[ic] += (area / 3.0) * sum_omega * (v_psi[ib] - v_psi[ia]);
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/* Your implementation goes here */

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

	// find psi_+1

	// find omega_+1

	t += dt;
}
