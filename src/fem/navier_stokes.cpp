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
	MnudtS = std::move(CSRMatrix(P,0.0));
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);

	build_P1_stiffness_matrix(m, MnudtS);
	blas_scal(nu*dt,MnudtS.data.data,MnudtS.data.size);
	build_P1_mass_matrix(m, MnudtS);

	S_chol = std::move(CholeskySolver(S));
	MnudtS_chol = std::move(CholeskySolver(S));
	
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

void NavierStokesSolver::compute_transport(double *T)
{
	// computes T(omega,psi)
	memset(T, 0, N * sizeof(double));
	for (size_t t=0;t<m.index_count();t+=3){
		uint32_t Ai = m.indices[t];
		uint32_t Bi = m.indices[t+1];
		uint32_t Ci = m.indices[t+2];
		uint32_t points[3] = {Ai,Bi,Ci};
		double sum = omega[Ai] + omega[Bi] + omega[Ci];
		for (int32_t k =0;k<3;k++)
			T[points[k]]+= sum * (psi[points[(3+k+1)%3]] - psi[points[(3+k-1)%3]]) / 6;
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

	double rel_error;
	iter = conjugate_gradient_solve(
			S,
			Momega.data,
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
	compute_stream_function(); // Momega computed here
	compute_transport(p.data);
	if (dt != this->dt || nu != this->nu){
		this->dt = dt;
		this->nu = nu;
		blas_copy(S.data.data,MnudtS.data.data,MnudtS.data.size);
		blas_axpby(1.0,M.data.data,dt*nu,MnudtS.data.data,MnudtS.data.size);
		MnudtS_chol.update_same_pattern(MnudtS);
		LOG_MSG("Updated cholesky decomposition");
		// MnudtS.dump("A.txt");
		// MnudtS_chol.L.dump("L.txt");
		// MnudtS_chol.L_anti_transpose.dump("L_antitrans.txt");
	}

	blas_axpby(1,Momega.data,dt,p.data,N);
	// at this stage p = b = M * omega(t) + dt * T(Omega,Psi)(t)
	// p.dump("p.txt");
	MnudtS_chol.solve(p.data,omega.data);
	// omega.dump("omega.txt");
	set_zero_mean(omega.data);
	t += dt;
}
