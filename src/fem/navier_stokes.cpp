#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <iostream>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"
#include "conjugate_gradient.h"

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
	M.mvp(V, Ap.data);
	double sum_VM = blas_sum_in_place(Ap.data, N);
	for(size_t i = 0;  i < N; i ++){
		V[i] -= sum_VM/vol;
	}
}

void NavierStokesSolver::compute_transport(double *T)
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

size_t NavierStokesSolver::compute_stream_function()
{
	
	//before iterations
	double b2 = blas_dot(Momega.data, Momega.data , N);
	S.mvp(psi.data, r.data);
	blas_axpby(1, psi.data, -1, r.data, N);
	/* p_0 = r_0 */
	blas_copy(r.data, p.data, N);

	double r2 = blas_dot(r.data, r.data, N);
	double rel_error = sqrt(r2/b2);


	size_t iter = 0;

	while ((iter < iter_max) && (rel_error > tol)){
		S.mvp(p.data, Ap.data);

		float alpha = r2/blas_dot(p.data, Ap.data, N); //compute alpha
		blas_axpy(alpha, p.data, psi.data, N); //Compute new x
	 	blas_axpy(-alpha, Ap.data, r.data, N); // Compute new r
	 	double r2_new = blas_dot(r.data, r.data, N);
	
	 	double beta = r2_new/ r2;
	 	blas_axpby(1, p.data, beta , r.data, N); //compute new p
		r2 = r2_new;
		rel_error = sqrt(r2 / b2);
		iter++;

	}

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

	compute_transport(p.data); //we save it in p as we will overwrite it anyway later
	//we are putting the 'b'-term of the GC in T
	blas_axpby(1, Momega.data, dt, p.data, N);
	//before iterations
	double b2 = blas_dot(p.data, p.data , N);

	// for x_0 take omgea_n --> r_0 = Momega + v*t*S
	S.mvp(omega.data, r.data); //first compute Somgea
	blas_axpby(1, Momega.data, dt*nu, r.data, N); //then we add with Momega that we have
	blas_axpby(1, p.data, -1, r.data, N); //then compute r_0
	
	/* p_0 = r_0 */
	blas_copy(r.data, p.data, N);

	double r2 = blas_dot(r.data, r.data, N);
	double rel_error = sqrt(r2/b2);

	size_t iter = 0;

	while ((iter < iter_max) && (rel_error > tol)){
		//update Ap
		S.mvp(p.data, Ap.data); //first do Sp
		M.mvp(p.data, Momega.data); //same here, we use this as storage beacuse we do not need it
		blas_axpby(1, Momega.data, nu*dt, Ap.data, N);

		float alpha = r2/blas_dot(p.data, Ap.data, N); //compute alpha
		blas_axpy(alpha, p.data, omega.data, N); //Compute new x
	 	blas_axpy(-alpha, Ap.data, r.data, N); // Compute new r
	 	double r2_new = blas_dot(r.data, r.data, N);
	
	 	double beta = r2_new/ r2;
	 	blas_axpby(1, p.data, beta , r.data, N); //compute new p
		r2 = r2_new;
		rel_error = sqrt(r2 / b2);
		iter++;
	}
	

	set_zero_mean(omega.data);

	t += dt;
}
