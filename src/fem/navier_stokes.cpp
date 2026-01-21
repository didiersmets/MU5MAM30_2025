#include <assert.h>
#include <stddef.h>
#include <stdint.h>

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
	double sum_M = M.sum();
	double sum_VM;
	M.mvp(V, sum_VM);
	V -= sum_VM/sum_M;
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	size_t num_ind = m.index_count;

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
	//prepare everything

	double x[N];
	memset(X, 0, N * sizeof(double));
	double rel_error = 0;
	double 
	


	size_t iter = conjugate_gradient_solve(S, -Momega, x, r, p, Ap, rel_error, tol,  iter_max, inited );

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
	double T[N];

	compute_transport(T);

	double b[N] = Momega + t*T;

	double b2 = blas_dot(b, b, N);

	if (!inited) {
		/* r_0 = b - Ax_0 */
		A.mvp(x, r);
		blas_axpby(1, b, -1, r, N);
		/* p_0 = r_0 */
		blas_copy(r, p, N);
	}

	double r2 = blas_dot(r, r, N);
	*rel_error = sqrt(r2 / b2);

	int iter = 0;
	while ((iter < iter_max) && (*rel_error > tol)) {
		size_t N = A.rows;
	 
		A.mvp(p, Ap); // update Ap
	 	float alpha = r2 / blas_dot(p, Ap, N); //Compute alpha
	 	blas_axpy(alpha, p, x, N); //Compute new x
	 	blas_axpy(-alpha, Ap, r, N); // Compute new r
	 	double r2_new = blas_dot(r, r, N);
	
	 	double beta = r2_new/ r2;
	 	blas_axpby(1, p, beta , r, N); //compute new p
		*rel_error = sqrt(r2 / b2);
		iter++;
	}
	

	set_zero_mean(omega.data);

	t += dt;
}
