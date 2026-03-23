#include "poisson_non_homogene_Dirichlet.h"

#include "P1.h"
#include "array.h"
#include "mesh.h"
#include "tiny_blas.h"
#include "projected_conjugate_gradient.h"

#include <map>
#include <utility>

PoissonSolverDirichlet::PoissonSolverDirichlet(const Mesh &m, BoundaryFunc g)
    : m(m), N(m.vertex_count()), g_values(N), f(N), u(N, 0.0), r(N), p(N), Ap(N)
{
#if USE_FEM_MATRIX
    build_P1_mass_matrix(m, M);
    build_P1_stiffness_matrix(m, A);
#else
    build_P1_CSRPattern(m, P);
    build_P1_mass_matrix(m, P, M);
    build_P1_stiffness_matrix(m, P, A);
#endif
    vol       = M.sum();
    inited    = false;
    iterate   = 0;
    converged = false;

    // ---- identify boundary vertices ----------------------------------------
    is_boundary.assign(N, false);
    std::map<std::pair<uint32_t,uint32_t>, int> edge_count;
    for (size_t t = 0; t < m.triangle_count(); t++) {
        uint32_t a = m.indices[3*t+0];
        uint32_t b = m.indices[3*t+1];
        uint32_t c = m.indices[3*t+2];
        uint32_t tri[3][2] = {{a,b},{b,c},{c,a}};
        for (auto &e : tri) {
            uint32_t lo = e[0] < e[1] ? e[0] : e[1];
            uint32_t hi = e[0] < e[1] ? e[1] : e[0];
            edge_count[{lo, hi}]++;
        }
    }
    has_boundary = false;
    for (auto &kv : edge_count) {
        if (kv.second == 1) {
            is_boundary[kv.first.first]  = true;
            is_boundary[kv.first.second] = true;
            has_boundary = true;
        }
    }

    // ---- evaluate g on boundary, zero elsewhere ----------------------------
    for (size_t i = 0; i < N; i++) {
        if (is_boundary[i]) {
            double x = m.positions[i][0];
            double y = m.positions[i][1];
            double z = m.positions[i][2];
            g_values.data[i] = g(x, y, z);
        } else {
            g_values.data[i] = 0.0;
        }
    }

    // ---- initialise solution to g_h ----------------------------------------
    for (size_t i = 0; i < N; i++)
        u[i] = g_values.data[i];
}

// Only used when there is no boundary (no-boundary fallback, same as poisson.cpp)
void PoissonSolverDirichlet::set_zero_mean(double *V)
{
    if (has_boundary) return;
    M.mvp(V, Ap.data);
    double s = blas_sum_in_place(Ap.data, N);
    for (size_t i = 0; i < N; ++i)
        V[i] -= s / vol;
}

void PoissonSolverDirichlet::clear_solution()
{
    // Reset solution to boundary values
    for (size_t i = 0; i < N; i++)
        u[i] = g_values.data[i];
    init_cg();
    iterate   = 0;
    converged = false;
}

void PoissonSolverDirichlet::init_cg()
{
    double *F  = f.data;
    double *U  = u.data;
    double *R  = r.data;
    double *Pv = p.data;
    double *AP = Ap.data;

    // No boundary: enforce zero mean for existence & uniqueness
    if (!has_boundary) {
        set_zero_mean(F);
        set_zero_mean(U);
    }

    // R = M*f
    M.mvp(F, R);

    // R = M*f - A*u,  with u = g_h on boundary => residual already encodes the lift
    A.mvp(U, AP);
    blas_axpy(-1.0, AP, R, N);

    // Zero out boundary rows of the residual:
    // boundary DOFs are already fixed at g, so they must not drive CG updates
    if (has_boundary) {
        for (size_t i = 0; i < N; i++)
            if (is_boundary[i]) R[i] = 0.0;
    }

    // b2 = ||R_0||^2  (reference norm for rel_error)
    b2 = blas_dot(R, R, N);

    blas_copy(R, Pv, N);
    r2        = blas_dot(R, R, N);
    rel_error = (b2 > 0) ? sqrt(r2 / b2) : 0.0;

    inited = true;
}

void PoissonSolverDirichlet::do_iterate(size_t max_iter, double tol)
{
    if (!inited)
        init_cg();

    while (max_iter-- && rel_error > tol) {
        // PCG: boundary vertices are projected out each iteration
        r2 = projected_cg_iterate_once(
            A, u.data, r.data, p.data, Ap.data, r2,
            has_boundary, has_boundary ? &is_boundary : nullptr
        );
        rel_error = (b2 > 0) ? sqrt(r2 / b2) : 0.0;
        iterate++;
    }

    if (rel_error <= tol)
        converged = true;
}