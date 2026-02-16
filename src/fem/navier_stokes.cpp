#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "navier_stokes.h"

#include "P1.h"
#include "P2.h" // Nuevo
#include "fem_matrix.h" // Nuevo
#include "tiny_blas.h"
#include "fem_type.h"


NavierStokesSolver::NavierStokesSolver(const Mesh &m, FEMType fem_type)
    : m(m), N(0), omega(), Momega(), psi(), r(), p(), Ap()
{
    switch (fem_type) {
    case FEMType::P1:
        build_P1_CSRPattern(m, P);
        build_P1_mass_matrix(m, P, M);
        build_P1_stiffness_matrix(m, P, S);
        break;

    case FEMType::P2:
        build_P2_CSRPattern(m, P);
        build_P2_mass_matrix(m, P, M);
        build_P2_stiffness_matrix(m, P, S);
        break;
    }

    N = P.rows;

    omega.resize(N);
    Momega.resize(N);
    psi.resize(N);
    r.resize(N);
    p.resize(N);
    Ap.resize(N);

    vol = M.sum();
    inited = false;
    t = 0;
}


void NavierStokesSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i) {
		V[i] -= s / vol;
	}
}

void NavierStokesSolver::compute_transport(double *T)
{
    memset(T, 0, N * sizeof(double));

    size_t vtx_count = m.vertex_count();

    // Transporte solo en vértices
    for (size_t t = 0; t < m.triangle_count(); t++) {
        uint32_t a = m.indices[3 * t + 0];
        uint32_t b = m.indices[3 * t + 1];
        uint32_t c = m.indices[3 * t + 2];

        assert(a < vtx_count && b < vtx_count && c < vtx_count);

        double sum = omega[a] + omega[b] + omega[c];
        T[a] += sum * (psi[b] - psi[c]);
        T[b] += sum * (psi[c] - psi[a]);
        T[c] += sum * (psi[a] - psi[b]);
    }

    // Escalado solo en vértices
    for (size_t v = 0; v < vtx_count; v++) {
        T[v] *= 1.0 / 6.0;
    }

    // Los dofs de arista (si los hay) se quedan en 0
    for (size_t v = vtx_count; v < N; ++v) {
        T[v] = 0.0;
    }
}



size_t NavierStokesSolver::compute_stream_function()
{
	double b2, r2, rel_error;
	size_t iter;

	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;
	double *Psi = psi.data;

	M.mvp(Om, MOm);

	/* Compute rhs norm2 */
	b2 = blas_dot(MOm, MOm, N);

	/* Form initial R and P */
	S.mvp(Psi, R);
	blas_axpby(1, MOm, -1, R, N);
	blas_copy(R, P, N);
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	/* Iterate until convergence */
	iter = 0;
	while ((rel_error > tol) && (iter++ < iter_max)) {

		/* Compute AP */
		S.mvp(P, AP);

		/* Update Psi */
		double alpha = r2 / blas_dot(P, AP, N);
		blas_axpy(alpha, P, Psi, N);

		/* Update R */
		blas_axpy(-alpha, AP, R, N);

		/* Update r2 and P */
		double beta = 1.0 / r2;
		r2 = blas_dot(R, R, N);
		rel_error = sqrt(r2 / b2);
		beta *= r2;
		blas_axpby(1, R, beta, P, N);
	}

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	double b2, r2, rel_error;

	size_t iter1, iter2;

	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;

	iter1 = compute_stream_function();

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Form rhs, saved in P */
	compute_transport(P);
	blas_axpby(1, MOm, dt, P, N);
	b2 = blas_dot(P, P, N);

	/* Form initial R and P */
	S.mvp(Om, R);
	blas_axpby(1, MOm, dt * nu, R, N);
	blas_axpby(1, P, -1, R, N);
	blas_copy(R, P, N);
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	/* Iterate until convergence (and at least once) */
	iter2 = 0;
	do {

		/* Compute AP (invalidates Mom) */
		S.mvp(P, AP);
		M.mvp(P, MOm); /* MOm used as temp storage */
		blas_axpby(1, MOm, dt * nu, AP, N);

		/* Update Om */
		double alpha = r2 / blas_dot(P, AP, N);
		blas_axpy(alpha, P, Om, N);

		/* Update R */
		blas_axpy(-alpha, AP, R, N);

		/* Update r2 and P */
		double beta = 1.0 / r2;
		r2 = blas_dot(R, R, N);
		rel_error = sqrt(r2 / b2);
		beta *= r2;
		blas_axpby(1, R, beta, P, N);

		/* Update MOm */
		M.mvp(Om, MOm);

		iter2++;
	} while ((rel_error > tol) && (iter2 <= iter_max));

	set_zero_mean(omega.data);

	t += dt;

	(void)iter1;
	//printf("Iter 1 : %zu, Iter2 : %zu\n", iter1, iter2);
}
