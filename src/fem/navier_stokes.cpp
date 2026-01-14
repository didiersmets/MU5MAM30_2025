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
	/* Your implementation goes here*/
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for(size_t i = 0; i < N; ++i){
		V[i] -= s / vol;
	}	
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */
	for (size_t i = 0; i < m.triangle_count();++i){
		//Get traingle vertex indices
		uint32_t i0 = m.indices[3*i];
		uint32_t i1 = m.indices[3*i+1];
		uint32_t i2 = m.indices[3*i+2];
		//Get triangle vertex positions
		Vec3 p0 = m.positions[i0];
		Vec3 p1 = m.positions[i1];
		Vec3 p2 = m.positions[i2];
		//Compute the area of the triangle
		Vec3 e1 = p1 - p0;
		Vec3 e2 = p2 - p0;
		Vec3 cross = Vec3(e1.y * e2.z - e1.z * e2.y,
				  e1.z * e2.x - e1.x * e2.z,
				  e1.x * e2.y - e1.y * e2.x);
		double area = 0.5 * sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
		//Compute gradients of basis functions
		Vec3 grad_phi0 ;
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

	t += dt;
}
