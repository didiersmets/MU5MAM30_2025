#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"
#include "conjugate_gradient.h"
#include "logging.h"

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
	M = std::move(CSRMatrix(P,0.0));
	S = std::move(CSRMatrix(P,0.0));
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
	LOG_MSG("M :");
	// M.print();
	LOG_MSG("S :");
	// S.print();
#endif
	vol = M.sum();
	inited = true;
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
		double sum = omega[Ai] + omega[Bi] + omega[Ci];
		for (uint32_t k =0;k<3;k++)
			T[points[k]]+= dt * sum * (psi[points[(k+1)%3]] - psi[points[(k-1)%3]]) / 0.6;
		}
}

size_t NavierStokesSolver::compute_stream_function()
{
	// solves for psi in : S psi = -M omega
	size_t iter = 0;
	double omega_norm = blas_dot(omega.data,omega.data,N);
	LOG_MSG("compute_stream_function : ||omega|| = %lf",omega_norm);
	M.mvp(omega.data,Momega.data);
	double Momega_norm = blas_dot(Momega.data,Momega.data,N);
	LOG_MSG("compute_stream_function : ||Momega|| = %lf",Momega_norm);

	TArray<double> b (N,0.0); // TODO : optmization possible
	blas_axpy(-1.0,Momega.data,b.data,N);

	double rel_error;
	iter = conjugate_gradient_solve(
			S,
			b.data,
			psi.data,
			r.data,
			p.data,
			Ap.data,
			&rel_error,
			10e-6,
			500,
			true);
	LOG_MSG("compute stream function : %d iter; %lf rel_error",iter,rel_error);
	return iter;
}
double NavierStokesSolver::cg_iterate_once(
		double dt, double nu, double *__restrict x,
		double *__restrict r, double *__restrict p,
		double *__restrict Ap, double r2){
	// Ap = (M + nu dt S)p
	S.mvp(p,Ap);
	blas_scal(dt * nu,Ap,N);
	M.add_mvp(p,Ap);
	double alpha = r2/blas_dot(p,Ap,N); //alpha = r2/p2_A
	blas_axpy(alpha,p,x,N); // x = x + alpha*p
	blas_axpy(-alpha,Ap,r,N); // r = r -alpha*A*p
	double new_r2 = blas_dot(r,r,N); //r2_{n+1}
	double beta = new_r2/r2; // beta = r2_{n+1}/r2_n
	blas_axpby(1,r,beta,p,N); // p = r + beta p
	return new_r2;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/
	compute_stream_function(); // Momega computed here
	// compute_transport(Momega.data,dt); // TODO : possible optimization
	TArray<double> T(N,0.0);
	compute_transport(T.data,dt);

	TArray<double> b (N,0.0);
	blas_copy(Momega.data,b.data,N);
	blas_axpy(1.0,T.data,b.data,N);
	// at this stage b = M * omega(t) + dt * T(Omega,Psi)(t)
	// initialization
	// r0 = b - Ax0 = b - (M + dt * nu * S) Omega
	S.mvp(omega.data,r.data);
	blas_axpby(1.0,Momega.data,dt*nu,r.data,N);
	blas_axpby(1.0,b.data,-1.0,r.data,N);
	blas_copy(r.data,p.data,N);
	// code bellow copied from conjugate_gradient.cpp
	double b2 = blas_dot(b.data,b.data,N);
	LOG_MSG("time step : b2 = %lf",b2);
	double r2 = blas_dot(r.data,r.data,N);
	double rel_error = sqrt(r2/b2);
	size_t iter = 0;
	size_t max_iter = 500;
	TArray<double> x (N,0.0); // TODO : possible optimization solve directly in omega
	while(iter<max_iter && rel_error>tol){
		r2 = cg_iterate_once(dt, nu,x.data, r.data, p.data, Ap.data, r2);
		rel_error = r2/b2;
		iter++;
	}
	blas_copy(x.data,omega.data,N);

	LOG_MSG("time step : %d iterations; %f relative error",iter,rel_error);
	set_zero_mean(omega.data);
	t += dt;
}
