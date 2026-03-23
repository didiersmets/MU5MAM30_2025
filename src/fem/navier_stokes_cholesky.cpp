#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes_cholesky.h"

#include "P1.h"
#include "tiny_blas.h"

#include "conjugate_gradient.h"
#include "cholesky.h"
#include "vec3.h"

NavierStokesSolver::NavierStokesSolver(const Mesh &m, double dt, double nu)
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

	// ---- Factorize S (with diagonal shift to handle kernel) ----------------
	CSRMatrix S_shift;
	S_shift.symmetric = true;
	S_shift.rows      = N;
	S_shift.cols      = N;
	S_shift.nnz       = P.col.size;
	S_shift.row_start = P.row_start.data;
	S_shift.col       = P.col.data;
	S_shift.data.resize(P.col.size);
	for (size_t k = 0; k < P.col.size; ++k)
		S_shift.data.data[k] = S.data.data[k];
	double diag_sum = 0.0;
	for (uint32_t i = 0; i < (uint32_t)N; ++i)
		diag_sum += S.data.data[P.row_start[i+1] - 1];
	double shift = 1e-6 * diag_sum / N;
	for (uint32_t i = 0; i < (uint32_t)N; ++i)
		S_shift.data.data[P.row_start[i+1] - 1] += shift;
	chol_S.factorize(S_shift);

	// ---- Factorize A = M + nu*dt*S -----------------------------------------
	A_dt.symmetric  = true;
	A_dt.rows       = N;
	A_dt.cols       = N;
	A_dt.nnz        = P.col.size;
	A_dt.row_start  = P.row_start.data;
	A_dt.col        = P.col.data;
	A_dt.data.resize(P.col.size);
	for (size_t k = 0; k < P.col.size; ++k)
		A_dt.data.data[k] = M.data.data[k] + nu * dt * S.data.data[k];
	chol_A.factorize(A_dt);
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = 0.0;
	for (size_t i = 0; i < N; i++) s += Ap.data[i];
	double c = s / vol;
	for (size_t i = 0; i < N; i++) V[i] -= c;
}

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
	set_zero_mean(Momega.data);
	chol_S.solve(Momega.data, psi.data);
	set_zero_mean(psi.data);
	return 1;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();

	double *AP = Ap.data;
	double *Om = omega.data;

	compute_transport(AP);
	M.mvp(Om, Momega.data);
	blas_axpy(dt, AP, Momega.data, N);

	chol_A.solve(Momega.data, Om);

	set_zero_mean(omega.data);
	t += dt;
}