#include "poisson.h"

#include "P1.h"
#include "array.h"
#include "conjugate_gradient.h"
#if USE_FEM_MATRIX
#include "fem_matrix.h"
#else
#include "sparse_matrix.h"
#endif
#include "mesh.h"
#include "tiny_blas.h"

PoissonSolver::PoissonSolver(Mesh &m, const bool &use_fem_P2)
	: m(m), f(0), u(0, 0.0), r(0), p(0), Ap(0), use_fem_P2(use_fem_P2)
{
#if USE_FEM_MATRIX
	if (use_fem_P2)
		throw std::runtime_error("P2 Lagrange elements are not available for FEM matrix framework");
	else
	{
		build_P1_mass_matrix(m, M);
		build_P1_stiffness_matrix(m, A);
	}
#else
	if (!use_fem_P2)
	{
		if (!m.is_periodic)
		{
			build_P1_CSRPattern(m, P);
			build_P1_mass_matrix(m, P, M);
			build_P1_stiffness_matrix(m, P, A);
		}
		else
		{
			build_P1_CSRPattern_per(m, P);
			build_P1_mass_matrix_per(m, P, M);
			build_P1_stiffness_matrix_per(m, P, A);
		}
	}
	else
	{
		m.build_edges();
		build_P2_CSRPattern(m, P);
		build_P2_mass_matrix(m, P, M);
		build_P2_stiffness_matrix(m, P, A);
	}
#endif

	vol = M.sum();
	inited = false;
	iterate = 0;
	converged = false;

	N = M.cols;
	f.resize(N);
	u.resize(N);
	r.resize(N);
	p.resize(N);
	Ap.resize(N);
}

void PoissonSolver::clear_solution()
{
	for (size_t i = 0; i < N; i++)
	{
		u[i] = 0;
	}
	init_cg();
	iterate = 0;
	converged = false;
}

void PoissonSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i)
	{
		V[i] -= s / vol;
	}
}

void PoissonSolver::init_cg()
{
	double *F = f.data;
	double *U = u.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	/* Fix up F and U for zero mean */
	set_zero_mean(F);
	set_zero_mean(U);

	/* r_0 = Mf - Au_0 */
	/* Compute also b2 = ||Mf||^2  */
	M.mvp(F, R);
	b2 = blas_dot(R, R, N);
	A.mvp(U, AP); /* AP used as temp storage */
	blas_axpy(-1, AP, R, N);
	/* p_0 = r_0 */
	blas_copy(R, P, N);
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	inited = true;
}

void PoissonSolver::do_iterate(size_t max_iter, double tol)
{
	if (!inited)
	{
		init_cg();
	}

	double *U = u.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	while (max_iter-- && rel_error > tol)
	{
		r2 = cg_iterate_once(A, U, R, P, AP, r2);
		iterate++;
		rel_error = sqrt(r2 / b2);
	}

	if (rel_error <= tol)
		converged = true;
}
