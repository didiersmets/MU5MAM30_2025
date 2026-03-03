#include "poisson.h"

#if USE_P2
  #include "P2.h"
#else
  #include "P1.h"
#endif
#include "array.h"
#include "conjugate_gradient.h"
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif
#include "mesh.h"
#include "tiny_blas.h"

static size_t inline dof_count(const Mesh &m) {
#if USE_P2
  return  m.vertex_count() + m.index_count() / 2;
#else
  return  m.vertex_count();
 #endif
}

static inline void unpack(uint64_t p, uint32_t &a, uint32_t &b) noexcept {
    a = uint32_t(p >> 32);
    b = uint32_t(p & 0xFFFFFFFFu);
}

PoissonSolver::PoissonSolver(const Mesh &m)
  : m(m), N(dof_count(m)), f(N), u(N, 0.0), r(N), p(N), Ap(N)
{
#if USE_FEM_MATRIX
  build_P1_mass_matrix(m, M);
  build_P1_stiffness_matrix(m, A);
#elif USE_P2
  std::unordered_map<uint64_t, uint32_t> edge2dof;
  build_P2_CSRPattern(m, P, edge2dof);
  build_P2_mass_matrix(m, P, M, edge2dof);
  build_P2_stiffness_matrix(m, P, A, edge2dof);

  uint32_t n_vtx = m.vertex_count();
  uint32_t n_edges =  m.index_count() / 2;
  e2vtx.resize(2 * n_edges);
  for (const auto &elt : edge2dof) {
    uint32_t v0, v1;
    unpack(elt.first, v0, v1);
    uint32_t e = elt.second - n_vtx;
    e2vtx[2*e] = v0;
    e2vtx[2*e+1] = v1;
  }
#else
  N = m.vertex_count();
  build_P1_CSRPattern(m, P);
  build_P1_mass_matrix(m, P, M);
  build_P1_stiffness_matrix(m, P, A);
#endif
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
