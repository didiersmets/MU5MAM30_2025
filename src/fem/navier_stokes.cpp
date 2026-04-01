#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

#include "conjugate_gradient.h"

/*********************************************************
 * Uses the following implementations from the solutions:
 * - set_zero_mean
 * - compute_transport, even though the calculated terms
 * 	 are a bit different
 * - time_step (own version corrected by AI)
 * 
**********************************************************/

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

/* Member function */
/* sum(M*V') = sum(M*V) - sum(M*ONES)*s/vol = s-vol*s/vol = s-s = 0 */
void NavierStokesSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i) {
		V[i] -= s / vol; /* V' = V-s/vol*ONES */
	}
}

/* Member function */
void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	double coriolis = 1e10;
	for (size_t t = 0; t < m.triangle_count(); t++) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		assert(a < N && b < N && c < N);
		double sum = omega[a] + omega[b] + omega[c];
		double sumPointsZ = m.positions[a].z + m.positions[b].z + m.positions[c].z;
		T[a] += (sum + coriolis*sumPointsZ) * (psi[b] - psi[c]);
		T[b] += (sum + coriolis*sumPointsZ) * (psi[c] - psi[a]);
		T[c] += (sum + coriolis*sumPointsZ) * (psi[a] - psi[b]);
	}

	for (size_t v = 0; v < N; v++) {
		T[v] *= 1.0 / 6;
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	double* minMOm = Momega.data;
	double* temp = minMOm;
	blas_axpby(0, temp, -1, minMOm, N);

	double rel_error_val = 1;
	double* rel_error = &rel_error_val;

	conjugate_gradient_solve(S, minMOm, psi.data, r.data, p.data, Ap.data, rel_error, tol, iter_max, true);

	temp = minMOm;
	blas_axpby(0, temp, -1, minMOm, N);

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
    double *R = r.data;
    double *P = p.data;
	double *AP = Ap.data;
    double *Om = omega.data;
    double *MOm = Momega.data;

	compute_stream_function();

    compute_transport(P);

    M.mvp(Om, MOm);
    blas_axpby(1, MOm, dt, P, N);

    double b2 = blas_dot(P, P, N);

    S.mvp(Om, R);
    blas_axpby(1, MOm, dt * nu, R, N);
    blas_axpby(1, P, -1, R, N);

    blas_copy(R, P, N);

    double r2 = blas_dot(R, R, N);
    double rel_error = sqrt(r2 / b2);
    size_t iter = 0;

    do {
        S.mvp(P, AP);
        M.mvp(P, MOm);
        blas_axpby(1, MOm, dt * nu, AP, N);

        double alpha = r2 / blas_dot(P, AP, N);

        blas_axpy(alpha, P, Om, N);

        blas_axpy(-alpha, AP, R, N);

        double beta = 1.0 / r2;
        r2 = blas_dot(R, R, N);
        rel_error = sqrt(r2 / b2);
        beta *= r2;

        blas_axpby(1, R, beta, P, N);

        M.mvp(Om, MOm);
        iter++;

    } while ((rel_error > tol) && (iter <= iter_max));

    set_zero_mean(omega.data);
    t += dt;
}
