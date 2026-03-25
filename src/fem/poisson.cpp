#include "poisson.h"

#include "adjacency_P2.h"
#include "P1.h"
#include "array.h"
#include "conjugate_gradient.h"
#include "P2.h"
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif
#include "mesh.h"
#include "tiny_blas.h"
/*ajout de la prise en compte du degré*/
PoissonSolver::PoissonSolver(const Mesh &m, int degre)
    : m(m), N(m.vertex_count()), f(N), u(N, 0.0), r(N), p(N), Ap(N)
{
	if (degre == 1) {
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, A);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, A);
#endif
}
	else if (degre == 2) {
		EdgeAdjacency edge_adj(m);
		build_P2_CSRPattern(m, P, edge_adj);

		N = P.rows; // Mise à jour de la taille pour P2
		f.resize(N); u.resize(N); r.resize(N); p.resize(N); Ap.resize(N);
		for (size_t i = 0; i < N; ++i) u[i] = 0.0;

		// Assemblage des matrices pareil que P1
		build_P2_mass_matrix(m, P, M, edge_adj);
		build_P2_stiffness_matrix(m, P, A, edge_adj);
	}
	vol = M.sum();
	inited = false;
	iterate = 0;
	converged = false;
}

void PoissonSolver::clear_solution()
{
	for (size_t i = 0; i < N; i++) {
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
	for (size_t i = 0; i < N; ++i) {
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
	if (!inited) {
		init_cg();
	}

	double *U = u.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	while (max_iter-- && rel_error > tol) {
		r2 = cg_iterate_once(A, U, R, P, AP, r2);
		iterate++;
		rel_error = sqrt(r2 / b2);
	}

	if (rel_error <= tol)
		converged = true;
}
