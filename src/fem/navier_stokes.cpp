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
double NavierStokesSolver::cg_iterate_once(double dt, double nu, double *__restrict x,
		       double *__restrict r, double *__restrict p,
	       	double *__restrict Ap, double r2){
	size_t N = M.rows;
	// Ap = (M + nu dt S)p
	S.mvp(p,Ap);
	blas_scal(dt * nu);
	M.add_mvp(p,Ap);
	double alpha = r2/blas_dot(p,Ap,N); //alpha = r2/p2_A
	blas_axpy(alpha,p,omega,N); // x = x + alpha*p
	blas_axpy(-alpha,Ap,r,N); // r = r -alpha*A*p
	double new_r2 = blas_dot(r,r,N); //r2_{n+1}
	double beta = new_r2/r2; // beta = r2_{n+1}/r2_n
	blas_axpby(1,r,beta,p,N); // p = r + beta p
	return new_r2;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();
	compute_transport(Momega,dt);
	// at this stage Momega = b = M * omega(t) + dt * T(Omega,Psi)(t)
	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/
	// initialization
	// r0 = b - Ax0 = Momega - (M + dt * nu * S) Omega
	S.mvp(omega.data,r.data);
	blas_scal(dt * nu);
	M.add_mvp(p.data,r.data);
	blas_axpby(1,Momega,-1,r.data,N);
	blas_copy(r.data,p.data,N);
	// code bellow copied from conjugate_gradient.cpp
	double b2 = blas_dot(Momega,Momega,N);
	double r2 = blas_dot(r.data,r.data,N);
	*rel_error = r2/b2;
	size_t iter = 0;
	while(iter<max_iter && *rel_error>tol){
		r2 = cg_iterate_once(dt, nu,omega.data, r.data, p.data, Ap.data, r2);
		*rel_error = r2/b2;
		iter++;
	}
	set_zero_mean(omega.data);
	t += dt;
}
