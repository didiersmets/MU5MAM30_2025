#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

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
	/* Your implementation goes here*/
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i)
	{
		V[i] -= s / vol;
	}
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	// /* Your implementation goes here */
	// for (size_t i = 0; i < m.triangle_count(); ++i)
	// {
	// 	// Get traingle vertex indices
	// 	uint32_t i0 = m.indices[3 * i];
	// 	uint32_t i1 = m.indices[3 * i + 1];
	// 	uint32_t i2 = m.indices[3 * i + 2];
	// 	// Get triangle vertex positions
	// 	Vec3 p0 = m.positions[i0];
	// 	Vec3 p1 = m.positions[i1];
	// 	Vec3 p2 = m.positions[i2];
	// 	// Compute the area of the triangle
	// 	Vec3 e1 = p1 - p0;
	// 	Vec3 e2 = p2 - p0;
	// 	Vec3 cross = Vec3(e1.y * e2.z - e1.z * e2.y,
	// 					  e1.z * e2.x - e1.x * e2.z,
	// 					  e1.x * e2.y - e1.y * e2.x);
	// 	double area = 0.5 * sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
	// 	// Compute gradients of basis functions
	// 	Vec3 grad0 = Vec3(e2.y, -e2.x, 0);
	// 	Vec3 grad1 = Vec3(-e1.y, e1.x, 0);
	// 	Vec3 grad2 = Vec3(e1.y - e2.y, e2.x - e1.x, 0);
	// 	grad0 /= (2 * area);
	// 	grad1 /= (2 * area);
	// 	grad2 /= (2 * area);
	// 	// Commpute the transport terms
	// 	T[i0] = area / 3.0 * (omega[i0] * psi[i1] * dot(grad0, grad1) + omega[i0] * psi[i2] * dot(grad0, grad2) + omega[i1] * psi[i1] * dot(grad1, grad0) + omega[i1] * psi[i2] * dot(grad2, grad0) + omega[i2] * psi[i1] * dot(grad1, grad0) + omega[i2] * psi[i2] * dot(grad2, grad0));

	// 	T[i1] = area / 3.0 * (omega[i0] * psi[i0] * dot(grad0, grad0) + omega[i0] * psi[i2] * dot(grad2, grad1) + omega[i1] * psi[i0] * dot(grad0, grad1) + omega[i1] * psi[i2] * dot(grad2, grad1) + omega[i2] * psi[i0] * dot(grad0, grad1) + omega[i2] * psi[i1] * dot(grad1, grad1));

	// 	T[i2] = area / 3.0 * (omega[i0] * psi[i0] * dot(grad0, grad2) + omega[i0] * psi[i1] * dot(grad1, grad2) + omega[i1] * psi[i0] * dot(grad0, grad2) + omega[i1] * psi[i1] * dot(grad1,grad2) + omega[i2] * psi[i0] * dot(grad0,grad2) + omega[i2] * psi[i1] * dot(grad1,grad2));
	// }
	
	//C'était beaucoup trop lourd en terme de calcul pas besioin de calculer les gradients....

	for (size_t t = 0; t < m.triangle_count(); t++) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		assert(a < N && b < N && c < N);
		double sum = omega[a] + omega[b] + omega[c];
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
	double b2, r2, r2_temp, rel_error;
	size_t iter = 0;
	TArray<double> p(N);
	/* Your implementation goes here */
	M.mvp(omega.data, Momega.data);
	/* Compute rhs norm2 */
	b2 = blas_dot(Momega.data, Momega.data, N);
	/* Form initial R and P */
	S.mvp(psi.data, r.data);
	blas_axpby(1, Momega.data, -1, r.data, N);
	blas_copy(r.data, p.data, N);
	r2 = blas_dot(r.data, r.data, N);
	rel_error = sqrt(r2 / b2);

	/* Iterate until convergence */
	for(size_t i = 0; i < N; ++i)
	{
		S.mvp(p.data, Ap.data);
		double alpha = r2 / blas_dot(p.data, Ap.data, N);
		blas_axpy(alpha, p.data, psi.data, N); //Update psi
		blas_axpy(-alpha, Ap.data, r.data, N); //Residual

		r2_temp = blas_dot(r.data, r.data, N);
		double beta = r2_temp / r2;
		blas_axpby(1, r.data, beta, p.data, N);
		r2 = r2_temp;
		rel_error = sqrt(r2 / b2);

		if(rel_error < tol)
		{
			break;
		}
	}
	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	double b2, r2, rel_error;
	size_t iter = 0;
	compute_stream_function();
	
	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Your implementation goes here */
	
	// Compute transport term
	compute_transport(p.data);
	//Compute rhs = M * omega(t) + dt * T(Omega,Psi)(t)
	M.mvp(omega.data, Momega.data);
	blas_axpy(dt, p.data, Momega.data, N);
	b2 = blas_dot(Momega.data, Momega.data, N);

	// Form initial residual and search direction
	S.mvp(omega.data, r.data);
	blas_axpby(1, Momega.data, -1, r.data, N);
	blas_copy(r.data, p.data, N);
	r2 = blas_dot(r.data,r.data,N);
	rel_error = sqrt(r2 / b2);

	//Ierate until convergence
	for(size_t i = 0; i < iter_max; ++i){
		S.mvp(p.data, Ap.data);
		M.mvp(p.data, Ap.data);
		blas_axpy(nu * dt, Momega.data, Ap.data, N);

		double alpha = r2 / blas_dot(p.data, Ap.data, N);
		blas_axpy(alpha, p.data, omega.data, N); //Update omega
		blas_axpy(-alpha, Ap.data, r.data, N); //Update residual

		double r2_temp = blas_dot(r.data, r.data, N);
		double beta = r2_temp / r2;
		blas_axpby(1, r.data, beta, p.data, N);
		r2 = r2_temp;
		rel_error = sqrt(r2 / b2);

		if(rel_error < tol){
			break;
		}
	}

	set_zero_mean(omega.data);

	t += dt;
}
