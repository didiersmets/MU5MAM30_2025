#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "conjugate_gradient.h"
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
	/*double mean = 0.0;
	for (size_t i = 0; i < N; i++) {
		mean += V[i];
	}
	mean /= N;
	for (size_t i = 0; i < N; i++) {
		V[i] -= mean;
	}
	*/

	std::vector<double> ones(N, 1.0);
    std::vector<double> Mones(N), MV(N);

    // Mones = M * ones
    M.multiply(ones, Mones);
	M.multiply(V, MV);

    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        num += ones[i] * MV[i];  // = (1^T M omega)
        den += ones[i] * Mones[i];   // = (1^T M 1)
    }

    double mean = num / den;

    // omega <- omega - mean
    for (std::size_t i = 0; i < N; ++i)
        V[i] -= mean;

}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */
	for (size_t k = 0; k < N; k++) {
		for (size_t j = 0; j < N; j++) {
			for (size_t i = 0; i < N; i++) {
				T[k] += omega[j] * psi[i] /12.0;
			}
		}
		T[k] *= vol;
	}

}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/* Your implementation goes here */
	// Solve S * psi = M * omega
	M.multiply(omega.data, Momega.data);	
	// Conjugate gradient
	
	iter = conjugate_gradient_solve(S, Momega.data, psi.data, r.data, p.data, Ap.data,
				&rel_error, 1e-8, 1000, inited);
	
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

	// Left-hand side matrix: A = M + nu * dt * S
	FEMatrix A = M;
	for (size_t i = 0; i < A.nnz; i++) {
		A.data[i] += nu * dt * S.data[i];		

	}

	// Right-hand side vector: b = M * omega(t) + dt * T(omega, psi)(t)
	double *b = new double[N];
	M.multiply(omega.data, b);
	double *T = new double[N];
	compute_transport(T);
	blas_axpy(dt, T, b, N);


	// Solve the linear system A * omega(t+dt) = b
	size_t iter = conjugate_gradient_solve(A, b, omega.data, r.data, p.data, Ap.data,
				&rel_error, 1e-8, 1000, inited);

	
	set_zero_mean(omega.data);

	t += dt;



}
