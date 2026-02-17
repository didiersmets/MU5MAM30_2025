#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <iostream>

#include "navier_stokes_cholesky.h"

#include "P1.h"
#include "tiny_blas.h"
#include "symbolic.h"
#include "cholesky.h"

NavierStokesSolverCholesky::NavierStokesSolverCholesky(const Mesh &m, double dt, double nu)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
	, transport(N)

	
	
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);
#endif

	construct_etree( P, parent);
	construct_L_sparsity_pattern(P, L_pattern, parent);
	cholesky_factorization(S, L_pattern, L_S);
	
	CSRMatrix S_M;
	S_M.rows = S_M.cols = N;
	S_M.nnz = M.nnz; 
	S_M.row_start = M.row_start; 
	S_M.col = M.col; 
	S_M.data.resize(M.nnz);
	S_M.symmetric = true;
	
	for(size_t i = 0; i < M.nnz; i ++){		
			S_M.data[i] = M.data[i] + nu*dt*S.data[i];}

	
	cholesky_factorization(S_M, L_pattern, L_SM);

	vol = M.sum();
	t = 0;
	this->dt = dt;
	this->nu = nu;
}

void NavierStokesSolverCholesky::set_zero_mean(double *V)
{
	M.mvp(V, transport.data);
	double sum_VM = blas_sum_in_place(transport.data, N);
	for(size_t i = 0;  i < N; i ++){
		V[i] -= sum_VM/vol;
	}
}

void NavierStokesSolverCholesky::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	size_t num_ind = m.index_count();

	for(size_t i = 0; i < num_ind; i += 3){
		uint32_t A = m.indices[i];
		uint32_t B = m.indices[i+1];
		uint32_t C = m.indices[i+2];

		double omega_sum = (omega[A] + omega[B] + omega[C])/6; 

		T[A] += (psi[B] - psi[C])*omega_sum;
		T[B] += (psi[C] - psi[A])*omega_sum;
		T[C] += (psi[A] - psi[B])*omega_sum;
	}
	
}

void NavierStokesSolverCholesky::compute_stream_function()
{
	TArray<double> x(N);
	for(size_t i = 0; i < N; i ++){
		x[i] = Momega[i];
	}
	solve_cholesky(L_S, x.data, psi.data);
	
}


void NavierStokesSolverCholesky::time_step()
{
	M.mvp(omega.data, Momega.data);

	compute_stream_function();
	

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Your implementation goes here */

	compute_transport(transport.data); 
	blas_axpby(1, Momega.data, dt, transport.data, N);

	solve_cholesky(L_SM, transport.data, omega.data);

	set_zero_mean(omega.data);

	t += dt;
}
