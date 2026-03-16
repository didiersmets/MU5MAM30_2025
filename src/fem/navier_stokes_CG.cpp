#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <vector>
#include "matrix_sum.h"

#include "navier_stokes_CG.h"

#include "P1.h"
#include "tiny_blas.h"

#include "conjugate_gradient.h"

//this is the version en cours, i added "_CG" in title
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

// the original structures
// void NavierStokesSolver::set_zero_mean(double *V)
// {
// 	/* Your implementation goes here */
// }

// void NavierStokesSolver::compute_transport(double *T)
// {
// 	memset(T, 0, N * sizeof(double));

// 	/* Your implementation goes here */
// }

// size_t NavierStokesSolver::compute_stream_function()
// {
// 	size_t iter = 0;

// 	/* Your implementation goes here */

// 	return iter;
// }

// void NavierStokesSolver::time_step(double dt, double nu)
// {
// 	compute_stream_function();

// 	/**********************************************************************
// 	 * Solve the system :
// 	 *
// 	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
// 	 *
// 	 *********************************************************************/

// 	/* Your implementation goes here */

// 	set_zero_mean(omega.data);

// 	t += dt;
// }


//calculate the intrgral(use the calculation from 1 to base functions to M), then minus it to ensure existence & unicity
//existence: S is SPD if restricted to the orthogonal of its kernel of dimension 1 (corresp. to constant functions)
//unicity: ensure psi don't include c*1(vec) ant for omega in time_step
void NavierStokesSolver::set_zero_mean(double *V)
{
    //  vec1^T * M * V = ∫ V_h dA
    M.mvp(V, Ap.data);
    double s = 0.0;
    for (size_t i = 0; i < N; i++) s += Ap.data[i];

    // c = ∫V dA / vol
    double c = s / vol;

    // V <- V - c， so 1^T * M * V_new = 0
    for (size_t i = 0; i < N; i++) V[i] -= c;
}


//key point: cross multiply + elimination of the triangle's area
void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	for (size_t t = 0; t < m.triangle_count(); t++) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];

		double sum = omega[a] + omega[b] + omega[c];

		T[a] += sum * (psi[c] - psi[b]);
		T[b] += sum * (psi[a] - psi[c]);
		T[c] += sum * (psi[b] - psi[a]);
	}

	for (size_t v = 0; v < N; v++) T[v] *= (1.0 / 6.0);
}


size_t NavierStokesSolver::compute_stream_function()
{
    M.mvp(omega.data, Momega.data);
    for (size_t i = 0; i < N; i++)
        Momega.data[i] = -Momega.data[i];

    set_zero_mean(Momega.data); // existence: 1^T * Momega = 0, use the calculation from 1 to base functions to M

    double rel_error = 0.0;
    size_t iter = conjugate_gradient_solve(
        S,
        Momega.data,
        psi.data,
        r.data,
        p.data,
        Ap.data,
        &rel_error,
        tol,
        iter_max,
        false
    );

    set_zero_mean(psi.data); // unicity: 1^T * psi = 0

    return iter;
}



// 1. Compute psi^m  by solving  S * psi = -M * omega
// 2. b = M*omega + dt*T(omega, psi)
// 3. Run CG to solve  (M + nu*dt*S) * omega_new = b,
//    starting from the current omega as initial guess
// 4. Enforce zero mean on the new omega
void NavierStokesSolver::time_step(double dt, double nu)
{
    // Step 1: compute stream function
    compute_stream_function();

    // Step 2: assemble RHS b = M*omega + dt*T
    std::vector<double> T_vec(N, 0.0);
    std::vector<double> b_vec(N, 0.0);

    compute_transport(T_vec.data());                    // T = T(omega, psi)
    M.mvp(omega.data, b_vec.data());                    // b = M*omega
    blas_axpy(dt, T_vec.data(), b_vec.data(), N);       // b += dt*T

    // Step 3: CG solve (M + nu*dt*S) * omega_new = b
    MatrixSum A(M, S, nu * dt, N); //this is a function in matrix_sum.h to combine these elements
    double rel_error = 0.0;
    conjugate_gradient_solve(
        A,
        b_vec.data(),
        omega.data,
        r.data,
        p.data,
        Ap.data,
        &rel_error,
        tol,
        iter_max,
        false
    );

    set_zero_mean(omega.data); //unicity of omega_new
    t += dt;
}