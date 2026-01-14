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
		Vec3 grad0 = Vec3(e2.y, -e2.x,0);
		Vec3 grad1 = Vec3(-e1.y, e1.x, 0);
		Vec3 grad2 = Vec3(e1.y - e2.y, e2.x - e1.x, 0);
		grad0 /= (2 * area);
		grad1 /= (2 * area);
		grad2 /= (2 * area);
		//Compute vorticity at triangle centroid
		double omega0 = omega[i0];
		double omega1 = omega[i1];
		double omega2 = omega[i2];
		
		for(size_t j = 0; j < 3; ++j){
			uint32_t vi = m.indices[3*i + j];
			double psi_val = psi[vi];
			Vec3 grad_psi;
			if (j == 0) grad_psi = grad0;
			else if (j == 1) grad_psi = grad1;
			else grad_psi = grad2;
			//Velocity u = curl(psi)
			Vec3 u = Vec3(grad_psi.y, -grad_psi.x, 0);
			//Contribution to transport term
			T[vi] += omega0 * (u.x * grad0.x + u.y * grad0.y) * area / 3.0;
			T[vi] += omega1 * (u.x * grad1.x + u.y * grad1.y) * area / 3.0;
			T[vi] += omega2 * (u.x * grad2.x + u.y * grad2.y) * area / 3.0;
		}		
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
