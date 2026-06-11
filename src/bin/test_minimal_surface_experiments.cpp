/* Numerical experiments for the minimal surface solver.
 *
 * Usage:
 *   ./build/debug/test_minimal_surface_experiments --experiment N
 *   N = 1  convergence analysis (L2 / H1 errors, Scherk)
 *   N = 2  Newton convergence rate (warm-start perturbation)
 *   N = 3  sensitivity to initial guess
 *   N = 4  condition number along Newton path
 *   N = 5  condition number scaling with mesh refinement
 *   N = 6  effect of boundary data magnitude
 *   N = 7  energy decrease monitoring
 *   N = 8  two sources of ill-conditioning (Scherk family)
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <chrono>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>

#include "minimal_graph.h"
#include "P1.h"
#include "conjugate_gradient.h"
#include "tiny_blas.h"
#include "square.h"
#include "disk.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "array.h"
#include "vec2.h"
#include "vec3.h"
#include "mesh_bounds.h"


using Vec2d = TVec2<double>;

/* =========================================================================
 * Local geometry helpers
 * ========================================================================= */
static double tri_area_2d(const Vec2d &AB, const Vec2d &AC)
{
    return 0.5 * std::fabs(AB.x * AC.y - AB.y * AC.x);
}

static Vec2d grad_uh(const Vec2d &AB, const Vec2d &AC,
                     double u_a, double u_b, double u_c)
{
    double det = AB.x * AC.y - AB.y * AC.x;
    if (std::fabs(det) < 1e-30) return {0.0, 0.0};
    double du_b = u_b - u_a, du_c = u_c - u_a;
    return {( AC.y * du_b - AB.y * du_c) / det,
            (-AC.x * du_b + AB.x * du_c) / det};
}

static Vec2d xy2d(const Mesh &m, uint32_t i)
{
    return {static_cast<double>(m.positions[i].x),
            static_cast<double>(m.positions[i].y)};
}

/* =========================================================================
 * Error computation (L2 and H1 seminorm)
 * ========================================================================= */
static void compute_errors(const Mesh &m, const TArray<double> &u_h,
                            std::function<double(const Vec2d &)> u_exact,
                            std::function<Vec2d(const Vec2d &)> grad_exact,
                            double &L2err, double &H1err)
{
    double L2sq = 0.0, H1sq = 0.0;
    for (size_t t = 0; t < m.triangle_count(); ++t) {
        uint32_t a = m.indices[3*t], b = m.indices[3*t+1], c = m.indices[3*t+2];
        Vec2d A = xy2d(m,a), B = xy2d(m,b), C = xy2d(m,c);
        Vec2d AB = {B.x-A.x, B.y-A.y}, AC = {C.x-A.x, C.y-A.y};
        double area = tri_area_2d(AB, AC);
        double ea = u_exact(A)-u_h[a], eb = u_exact(B)-u_h[b], ec = u_exact(C)-u_h[c];
        L2sq += area * (ea*ea + eb*eb + ec*ec) / 3.0;
        Vec2d gu = grad_uh(AB, AC, u_h[a], u_h[b], u_h[c]);
        Vec2d cen = {(A.x+B.x+C.x)/3.0, (A.y+B.y+C.y)/3.0};
        Vec2d ge  = grad_exact(cen);
        double dx = ge.x-gu.x, dy = ge.y-gu.y;
        H1sq += area * (dx*dx + dy*dy);
    }
    L2err = std::sqrt(L2sq);
    H1err = std::sqrt(H1sq);
}

/* =========================================================================
 * Boundary mask helper
 * ========================================================================= */
static void build_boundary_mask(const Mesh &m, std::vector<bool> &is_bnd)
{
    is_bnd.assign(m.vertex_count(), false);
    for (size_t i = 0; i < m.boundary.size; ++i)
        is_bnd[m.boundary[i]] = true;
}

static void rescale_and_recenter_mesh(Mesh &mesh)
{
    Aabb bbox = compute_mesh_bounds(mesh);
    Vec3 model_center = (bbox.min + bbox.max) * 0.5f;
    Vec3 model_extent = (bbox.max - bbox.min);
    float model_size = max(model_extent);
    if (model_size == 0) {
        printf("Warning : Mesh is empty or reduced to a point.\n");
        model_size = 1;
    }
    for (size_t i = 0; i < mesh.vertex_count(); ++i) {
        mesh.positions[i] -= model_center;
        mesh.positions[i] /= (model_size / 2);
    }
}

/* =========================================================================
 * Condition number estimate (interior block of J_F)
 *
 * Strategy:
 *   1. Power iteration with boundary DOFs projected to 0 -> lambda_max
 *      of the interior block (v_B=0 => v^T A v = v_I^T A_II v_I).
 *   2. Inverse power iteration via CG -> lambda_min of interior block.
 *   3. kappa = lambda_max / lambda_min.
 * ========================================================================= */
static double estimate_condition_number(const CSRMatrix &A, const Mesh &m,
                                        size_t N_pow = 60)
{
    size_t N = m.vertex_count();
    std::vector<bool> is_bnd;
    build_boundary_mask(m, is_bnd);

    size_t N_int = 0;
    for (size_t i = 0; i < N; ++i) if (!is_bnd[i]) ++N_int;
    if (N_int == 0) return 1.0;

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    /* ---- lambda_max via projected power iteration ---- */
    TArray<double> v(N, 0.0), w(N, 0.0);
    for (size_t i = 0; i < N; ++i) if (!is_bnd[i]) v[i] = dist(rng);
    for (size_t i = 0; i < N; ++i) if ( is_bnd[i]) v[i] = 0.0;
    double vnm = std::sqrt(blas_dot(v.data, v.data, N));
    if (vnm > 0) for (size_t i = 0; i < N; ++i) v[i] /= vnm;

    double lambda_max = 1.0;
    for (size_t k = 0; k < N_pow; ++k) {
        A.mvp(v.data, w.data);
        for (size_t i = 0; i < N; ++i) if (is_bnd[i]) w[i] = 0.0;
        lambda_max = blas_dot(v.data, w.data, N);
        double wn = std::sqrt(blas_dot(w.data, w.data, N));
        if (wn < 1e-30) break;
        for (size_t i = 0; i < N; ++i) v[i] = w[i] / wn;
    }

    /* ---- lambda_min via inverse power iteration + CG ---- */
    for (size_t i = 0; i < N; ++i) v[i] = is_bnd[i] ? 0.0 : dist(rng);
    vnm = std::sqrt(blas_dot(v.data, v.data, N));
    if (vnm > 0) for (size_t i = 0; i < N; ++i) v[i] /= vnm;

    TArray<double> y(N, 0.0), r_cg(N, 0.0), p_cg(N, 0.0), Ap_cg(N, 0.0);

    double lambda_min = lambda_max;
    for (size_t inv_k = 0; inv_k < 5; ++inv_k) {
        memset(y.data, 0, N * sizeof(double));
        double cg_err = 0.0;
        conjugate_gradient_solve(A, v.data, y.data,
                                 r_cg.data, p_cg.data, Ap_cg.data,
                                 &cg_err, 1e-10, 5000, false);
        for (size_t i = 0; i < N; ++i) if (is_bnd[i]) y[i] = 0.0;
        double yn = std::sqrt(blas_dot(y.data, y.data, N));
        if (yn < 1e-30) break;
        A.mvp(v.data, w.data);
        for (size_t i = 0; i < N; ++i) if (is_bnd[i]) w[i] = 0.0;
        lambda_min = blas_dot(v.data, w.data, N);
        for (size_t i = 0; i < N; ++i) v[i] = y[i] / yn;
    }

    if (lambda_min < 1e-30) return 1e30;
    return lambda_max / lambda_min;
}

/* =========================================================================
 * Helper: initialise CSRMatrix from existing pattern
 * ========================================================================= */
static void init_csr_from_pattern(const CSRPattern &P, CSRMatrix &S)
{
    S.symmetric = true;
    S.rows = S.cols = P.rows;
    S.nnz  = P.nnz;
    S.row_start = P.row_start.data;
    S.col       = P.col.data;
    S.data.resize(S.nnz);
    for (size_t i = 0; i < S.nnz; ++i) S.data[i] = 0.0;
}

/* =========================================================================
 * Instrumented Newton iteration
 * ========================================================================= */
struct NewtonRecord {
    size_t iter;
    double residual_norm;
    double correction_norm;
    double area;
    double alpha;
};

static bool run_newton_instrumented(
    const Mesh &m,
    std::function<double(const Vec2d &)> boundary_func,
    TArray<double> &u_init,
    size_t max_iter, double tol,
    std::vector<NewtonRecord> &records,
    double min_alpha = 1.0)
{
    size_t N = m.vertex_count();
    std::vector<bool> is_bnd;
    build_boundary_mask(m, is_bnd);

    CSRPattern P;
    build_P1_CSRPattern(m, P);

    TArray<double> u(N), du(N,0.0), q(m.triangle_count());
    TArray<double> rhs(N,0.0), u_tmp(N,0.0);
    TArray<double> r_cg(N,0.0), p_cg(N,0.0), Ap_cg(N,0.0);
    CSRMatrix S;
    init_csr_from_pattern(P, S);

    memcpy(u.data, u_init.data, N*sizeof(double));

    constexpr double HUGE_DIAG = 1e30;
    constexpr double c_armijo  = 1e-4;
    constexpr double rho       = 0.5;
    const double     tolCG     = 1e-14;

    MinimalGraphSolver helper(m, boundary_func);
    double area = helper.compute_denominator(q, u);
    bool converged = false;

    for (size_t k = 0; k < max_iter; ++k) {
        build_P1_stiffness_matrix_NS(m, P, S, q.data, u.data, area);
        build_P1_rhs_NS(m, q.data, u.data, rhs);

        double res_sq = 0.0;
        for (size_t i = 0; i < N; ++i)
            if (!is_bnd[i]) res_sq += rhs[i]*rhs[i];
        double res_norm = std::sqrt(res_sq);

        for (size_t i = 0; i < m.boundary.size; ++i) {
            uint32_t bi = m.boundary[i];
            S(bi,bi) = HUGE_DIAG;
            rhs[bi]  = 0.0;
        }

        memset(du.data, 0, N*sizeof(double));
        double cg_err = 0.0;
        conjugate_gradient_solve(S, rhs.data, du.data,
                                 r_cg.data, p_cg.data, Ap_cg.data,
                                 &cg_err, tolCG, 10000, false);

        double corr_sq = blas_dot(du.data, du.data, N);

        double alpha = 1.0;
        bool flag = true;
        while (flag) {
            for (size_t i = 0; i < N; ++i) u_tmp[i] = u[i] + alpha*du[i];
            double e_tmp = helper.compute_denominator(q, u_tmp);
            if (e_tmp < area - c_armijo*alpha*corr_sq || alpha <= min_alpha)
                flag = false;
            else
                alpha *= rho;
        }
        for (size_t i = 0; i < N; ++i) u[i] += alpha*du[i];
        area = helper.compute_denominator(q, u);

        NewtonRecord rec{k, res_norm, std::sqrt(corr_sq), area, alpha};
        records.push_back(rec);

        if (rec.correction_norm < tol) { converged = true; break; }
        if (k >= 10) {
            double prev = records[records.size()-11].correction_norm;
            if (rec.correction_norm > prev * 0.9999) break;
        }
    }
    memcpy(u_init.data, u.data, N*sizeof(double));
    return converged;
}

/* =========================================================================
 * Standard initial guess (boundary interpolated, interior zero)
 * ========================================================================= */
static void make_standard_init(const Mesh &m,
                                std::function<double(const Vec2d &)> f,
                                TArray<double> &u)
{
    size_t N = m.vertex_count();
    u.resize(N);
    memset(u.data, 0, N*sizeof(double));
    for (size_t i = 0; i < m.boundary.size; ++i) {
        uint32_t bi = m.boundary[i];
        u[bi] = f(xy2d(m, bi));
    }
}

/* =========================================================================
 * CSV writer
 * ========================================================================= */
static void write_csv(const std::string &fname,
                      const std::vector<std::string> &header,
                      const std::vector<std::vector<double>> &rows)
{
    std::ofstream f(fname);
    for (size_t j = 0; j < header.size(); ++j) {
        f << header[j];
        if (j+1 < header.size()) f << ',';
    }
    f << '\n';
    for (auto &row : rows) {
        for (size_t j = 0; j < row.size(); ++j) {
            f << row[j];
            if (j+1 < row.size()) f << ',';
        }
        f << '\n';
    }
    printf("Wrote %s\n", fname.c_str());
}

/* =========================================================================
 * CG iteration counter (for experiment 8)
 * ========================================================================= */
static int count_cg_iters(const CSRMatrix &A,
                           const std::vector<double> &b,
                           double tol = 1e-6, int max_iter = 200000)
{
    size_t n = A.rows;
    std::vector<double> x(n,0.0), r(b), p(b), Ap(n,0.0);
    double rr0 = 0.0;
    for (size_t i = 0; i < n; ++i) rr0 += r[i]*r[i];
    if (rr0 < 1e-100) return 0;
    const double tol2 = tol*tol*rr0;
    double rr = rr0;
    for (int k = 0; k < max_iter; ++k) {
        A.mvp(p.data(), Ap.data());
        double pAp = 0.0;
        for (size_t i = 0; i < n; ++i) pAp += p[i]*Ap[i];
        if (std::fabs(pAp) < 1e-100) return k;
        double alpha = rr / pAp;
        double rr_new = 0.0;
        for (size_t i = 0; i < n; ++i) {
            x[i] += alpha*p[i];
            r[i] -= alpha*Ap[i];
            rr_new += r[i]*r[i];
        }
        if (rr_new < tol2) return k+1;
        double beta = rr_new / rr;
        for (size_t i = 0; i < n; ++i) p[i] = r[i] + beta*p[i];
        rr = rr_new;
    }
    return max_iter;
}

/* =========================================================================
 * Scherk surface helpers (shared by experiments 1, 2, 8)
 *
 * Standard Scherk on [-1,1]^2:
 *   scale = pi/4
 *   u(x,y)  = (1/scale) * ln(cos(scale*x) / cos(scale*y))
 *   grad u  = (-tan(scale*x), tan(scale*y))
 *   M       = tan(scale) = 1  (at corners)
 *
 * Scherk family for experiment 8:
 *   u_alpha(x,y) = (1/alpha) * ln(cos(alpha*x) / cos(alpha*y))
 *   M = tan(alpha),  alpha in (0, pi/2)
 * ========================================================================= */
static constexpr double SCHERK_SCALE = M_PI / 4.0;

static double u_scherk_std(const Vec2d &p)
{
    double xc = SCHERK_SCALE * p.x;
    double yc = SCHERK_SCALE * p.y;
    return (1.0 / SCHERK_SCALE) * std::log(std::cos(xc) / std::cos(yc));
}

static Vec2d grad_scherk_std(const Vec2d &p)
{
    return {-std::tan(SCHERK_SCALE * p.x),
             std::tan(SCHERK_SCALE * p.y)};
}

/* =========================================================================
 * EXPERIMENT 1 — Convergence analysis (L2 and H1 errors)
 *
 * Three test cases:
 *   1. Linear u = x+y  on [0,1]^2   (P1-exact: errors == 0)
 *   2. Flat   u = 0    on disk R=1  (P1-exact: errors == 0)
 *   3. Scherk surface  on [-1,1]^2  (genuine O(h^2)/O(h) rates)
 *
 * For test case 3 the solver domain is [-1,1]^2 (build_square_mesh +
 * rescale_and_recenter_mesh). The Scherk function u_scherk_std is a
 * valid minimal surface on that domain with M = tan(pi/4) = 1.
 * ========================================================================= */
static void experiment_1()
{
    printf("\n=== Experiment 1: Convergence Analysis ===\n");

    std::vector<std::string> header = {
        "test_case","level","N_subdiv","h","N_vtx",
        "L2_error","H1_error","rate_L2","rate_H1"
    };
    std::vector<std::vector<double>> rows;

    /* ------------------------------------------------------------------ */
    /* Test case 1: linear u = x+y on [0,1]^2                             */
    /* ------------------------------------------------------------------ */
    auto u_lin  = [](const Vec2d &p) { return p.x + p.y; };
    auto gu_lin = [](const Vec2d &)  -> Vec2d { return {1.0, 1.0}; };

    printf("\nTest case 1: linear u=x+y (P1-exact; errors at machine epsilon)\n");
    printf("%-6s %-8s %-10s %-12s %-12s\n","level","N","h","L2 error","H1 error");

    double prev_L2=0, prev_H1=0, prev_h=0;
    for (size_t li = 0; li < 6; ++li) {
        size_t N = 2u << li;   /* 2,4,8,16,32,64 */
        Mesh m; build_square_mesh(&m, N, 1.0, 1.0);
        double h = 1.0 / (double)N;

        TArray<double> u_sol(m.vertex_count());
        for (size_t i = 0; i < m.vertex_count(); ++i) u_sol[i] = u_lin(xy2d(m,i));
        std::vector<NewtonRecord> recs;
        run_newton_instrumented(m, u_lin, u_sol, 500, 1e-12, recs, 1.0);

        double L2e, H1e;
        compute_errors(m, u_sol, u_lin, gu_lin, L2e, H1e);
        printf("%-6zu %-8zu %-10.4f %-12.3e %-12.3e\n", li, N, h, L2e, H1e);
        rows.push_back({1,(double)li,(double)N,h,(double)m.vertex_count(),
                        L2e,H1e,0.0,0.0});
    }

    /* ------------------------------------------------------------------ */
    /* Test case 2: flat u=0 on disk R=1                                   */
    /* ------------------------------------------------------------------ */
    auto u_zero  = [](const Vec2d &) { return 0.0; };
    auto gu_zero = [](const Vec2d &) -> Vec2d { return {0.0, 0.0}; };

    printf("\nTest case 2: flat u=0, disk R=1 (errors identically 0)\n");
    printf("%-6s %-8s %-10s %-12s %-12s\n","level","N","h","L2 error","H1 error");

    for (size_t li = 0; li < 5; ++li) {
        size_t N = 2u << li;   /* 2,4,8,16,32 */
        Mesh m; build_disk_mesh(&m, N, 1.0);
        double h = 2.0 * M_PI / (6.0 * N);

        TArray<double> u_sol(m.vertex_count(), 0.0);
        std::vector<NewtonRecord> recs;
        run_newton_instrumented(m, u_zero, u_sol, 500, 1e-12, recs, 1.0);

        double L2e, H1e;
        compute_errors(m, u_sol, u_zero, gu_zero, L2e, H1e);
        printf("%-6zu %-8zu %-10.4f %-12.3e %-12.3e\n", li, N, h, L2e, H1e);
        rows.push_back({2,(double)li,(double)N,h,(double)m.vertex_count(),
                        L2e,H1e,0.0,0.0});
    }

    /* ------------------------------------------------------------------ */
    /* Test case 3: Scherk on [-1,1]^2                                     */
    /* ------------------------------------------------------------------ */
    printf("\nTest case 3: Scherk minimal surface on [-1,1]^2 (rate_L2~2, rate_H1~1)\n");
    printf("%-6s %-8s %-10s %-12s %-12s %-12s %-8s %-8s\n",
           "level","N","h","Newton res","L2 error","H1 error","rate_L2","rate_H1");

    prev_L2=0; prev_H1=0; prev_h=0;
    std::vector<size_t> levels_sch = {8,16,32,64,128};
    for (size_t li = 0; li < levels_sch.size(); ++li) {
        size_t N = levels_sch[li];
        Mesh m; build_square_mesh(&m, N, 1, 1);
        rescale_and_recenter_mesh(m);
        double h = 2.0 / (double)N;

        TArray<double> u_sol(m.vertex_count());
        make_standard_init(m, u_scherk_std, u_sol);

        std::vector<NewtonRecord> recs;
        bool conv = run_newton_instrumented(m, u_scherk_std, u_sol,
                                            500, 1e-12, recs, 0.1);

        double L2e, H1e;
        compute_errors(m, u_sol, u_scherk_std, grad_scherk_std, L2e, H1e);

        double newton_res = recs.empty() ? 0.0 : recs.back().correction_norm;
        double rL2=0, rH1=0;
        if (li>0 && prev_L2>1e-30 && L2e>1e-30)
            rL2 = std::log(L2e/prev_L2) / std::log(h/prev_h);
        if (li>0 && prev_H1>1e-30 && H1e>1e-30)
            rH1 = std::log(H1e/prev_H1) / std::log(h/prev_h);

        printf("%-6zu %-8zu %-10.4f %-12.3e %-12.3e %-12.3e %-8.2f %-8.2f  [%s %zu iter]\n",
               li, N, h, newton_res, L2e, H1e, rL2, rH1,
               conv?"conv":"div", recs.size());
        rows.push_back({3,(double)li,(double)N,h,(double)m.vertex_count(),
                        L2e,H1e,rL2,rH1});
        prev_L2=L2e; prev_H1=H1e; prev_h=h;
    }

    write_csv("experiment_1_convergence.csv", header, rows);
}

/* =========================================================================
 * EXPERIMENT 2 — Newton convergence rate (warm-start perturbation)
 *
 * Two sets of runs written to the same CSV:
 *   init = "cold"      standard init (interior=0), O(1) from U*
 *   init = "warm_Xe"   U* + eps*xi, eps in {1e-1, 1e-2, 1e-4}
 *
 * Diagnostics:
 *   lin_ratio   |F^{k+1}|/|F^k|          -> 0 for quadratic
 *   quad_ratio  |F^{k+1}|/|F^k|^2        -> O(1) for quadratic
 *   loglog_slope log|F^{k+1}|/log|F^k|   -> 2 for quadratic
 * ========================================================================= */
static void experiment_2()
{
    printf("\n=== Experiment 2: Newton Convergence Rate ===\n");

    Mesh m; build_square_mesh(&m, 16, 1, 1);
    rescale_and_recenter_mesh(m);
    const size_t Nvtx = m.vertex_count();

    /* Scherk on [-1,1]^2 */
    auto u_scherk = [](const Vec2d &p) { return u_scherk_std(p); };

    /* Interior index list */
    std::vector<bool> is_bnd(Nvtx, false);
    for (size_t i = 0; i < m.boundary.size; ++i) is_bnd[m.boundary[i]] = true;
    std::vector<size_t> interior;
    for (size_t i = 0; i < Nvtx; ++i) if (!is_bnd[i]) interior.push_back(i);
    const size_t Nint = interior.size();

    /* Phase 1: converge to U* */
    TArray<double> u_star(Nvtx);
    make_standard_init(m, u_scherk, u_star);
    {
        std::vector<NewtonRecord> recs;
        bool ok = run_newton_instrumented(m, u_scherk, u_star, 500, 1e-13, recs, 1.0);
        printf("Phase 1 (converge to U*): %s in %zu iters, final |F|=%.3e\n",
               ok?"OK":"FAIL", recs.size(),
               recs.empty()?0.0:recs.back().residual_norm);
    }

    /* Diagnostics helper */
    std::vector<std::string> row_labels;
    std::vector<std::vector<double>> all_rows;

    auto record_run = [&](const std::string &label,
                          const std::vector<NewtonRecord> &recs) {
        for (size_t k = 0; k < recs.size(); ++k) {
            double lin=0, quad=0, slope=0;
            if (k>=1 && recs[k-1].residual_norm > 1e-16) {
                lin  = recs[k].residual_norm / recs[k-1].residual_norm;
                quad = recs[k].residual_norm /
                       (recs[k-1].residual_norm * recs[k-1].residual_norm);
            }
            if (k>=2 && recs[k].residual_norm>1e-16
                     && recs[k-1].residual_norm>1e-16
                     && recs[k-2].residual_norm>1e-16) {
                double lc=std::log(recs[k].residual_norm);
                double lp=std::log(recs[k-1].residual_norm);
                double lpp=std::log(recs[k-2].residual_norm);
                if (std::abs(lp-lpp)>1e-12) slope=(lc-lp)/(lp-lpp);
            }
            printf("%-5zu  %-12.4e  %-12.4e  %-10.4f  %-12.4e  %-12.4f\n",
                   recs[k].iter, recs[k].residual_norm, recs[k].correction_norm,
                   lin, quad, slope);
            row_labels.push_back(label);
            all_rows.push_back({(double)recs[k].iter, recs[k].residual_norm,
                                 recs[k].correction_norm, lin, quad, slope});
        }
    };

    /* Phase 2a: cold start */
    printf("\n--- Cold start ---\n");
    {
        TArray<double> u_cold(Nvtx);
        make_standard_init(m, u_scherk, u_cold);
        std::vector<NewtonRecord> recs;
        run_newton_instrumented(m, u_scherk, u_cold, 200, 1e-13, recs, 1.0);
        record_run("cold", recs);
    }

    /* Phase 2b: warm starts */
    std::mt19937 rng(42);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (double eps : {1e-1, 1e-2, 1e-4}) {
        char label[32];
        std::snprintf(label, sizeof(label), "warm_%.0e", eps);
        printf("\n--- Warm start eps=%.0e ---\n", eps);

        std::vector<double> xi(Nint);
        double norm2=0; for (double &v : xi) { v=nd(rng); norm2+=v*v; }
        double inv_n = 1.0/std::sqrt(norm2);

        TArray<double> u_pert(Nvtx);
        for (size_t i=0;i<Nvtx;++i) u_pert[i]=u_star[i];
        for (size_t k=0;k<Nint;++k) u_pert[interior[k]] += eps*xi[k]*inv_n;

        std::vector<NewtonRecord> recs;
        run_newton_instrumented(m, u_scherk, u_pert, 30, 1e-13, recs, 1.0);
        record_run(std::string(label), recs);
    }

    /* Write CSV */
    FILE *fp = fopen("experiment_2_newton_convergence.csv","w");
    fprintf(fp,"init,iter,residual_norm,correction_norm,"
               "lin_ratio,quad_ratio,loglog_slope\n");
    for (size_t i=0; i<all_rows.size(); ++i) {
        const auto &r=all_rows[i];
        fprintf(fp,"%s,%.0f,%.10e,%.10e,%.10e,%.10e,%.10f\n",
                row_labels[i].c_str(),r[0],r[1],r[2],r[3],r[4],r[5]);
    }
    fclose(fp);
    printf("\nWrote experiment_2_newton_convergence.csv\n");
}

/* =========================================================================
 * EXPERIMENT 3 — Sensitivity to initial guess
 * ========================================================================= */
static void experiment_3()
{
    printf("\n=== Experiment 3: Sensitivity to Initial Guess ===\n");

    Mesh m; build_square_mesh(&m, 10, 1.0, 1.0);
    size_t N = m.vertex_count();

    auto f = [](const Vec2d &p) {
        return std::sin(M_PI*p.x)*std::sinh(M_PI*p.y)/std::sinh(M_PI);
    };

    TArray<double> u_std(N);
    make_standard_init(m, f, u_std);

    std::vector<std::string> header = {"alpha_or_sigma","kind","converged","n_iter","final_residual"};
    std::vector<std::vector<double>> rows;

    printf("\nAlpha-interpolation:\n");
    for (double alpha : {0.0,0.2,0.4,0.6,0.8,1.0}) {
        TArray<double> u_init(N, 0.0);
        for (size_t i=0;i<N;++i) u_init[i]=alpha*u_std[i];
        std::vector<NewtonRecord> recs;
        bool conv = run_newton_instrumented(m,f,u_init,300,1e-12,recs,1.0);
        double fr = recs.empty()?-1.0:recs.back().residual_norm;
        printf("alpha=%4.2f  conv=%s  iter=%3zu  |F|=%.3e\n",
               alpha,conv?"yes":"no",recs.size(),fr);
        rows.push_back({alpha,0.0,conv?1.0:0.0,(double)recs.size(),fr});
    }

    printf("\nRandom perturbation:\n");
    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0,1.0);
    for (double sigma : {0.01,0.1,0.5,1.0,2.0,5.0}) {
        TArray<double> u_init(N);
        for (size_t i=0;i<N;++i) u_init[i]=u_std[i]+sigma*noise(rng);
        for (size_t i=0;i<m.boundary.size;++i) {
            uint32_t bi=m.boundary[i]; u_init[bi]=u_std[bi];
        }
        std::vector<NewtonRecord> recs;
        bool conv = run_newton_instrumented(m,f,u_init,300,1e-12,recs,1.0);
        double fr = recs.empty()?-1.0:recs.back().residual_norm;
        printf("sigma=%5.2f  conv=%s  iter=%3zu  |F|=%.3e\n",
               sigma,conv?"yes":"no",recs.size(),fr);
        rows.push_back({sigma,1.0,conv?1.0:0.0,(double)recs.size(),fr});
    }
    write_csv("experiment_3_initial_guess.csv", header, rows);
}

/* =========================================================================
 * EXPERIMENT 4 — Condition number along the Newton path
 * ========================================================================= */
static void experiment_4()
{
    printf("\n=== Experiment 4: Condition Number Along Newton Path ===\n");

    Mesh m; build_square_mesh(&m, 12, 1.0, 1.0);
    size_t N = m.vertex_count();

    auto f = [](const Vec2d &p) {
        return std::sin(M_PI*p.x)*std::sinh(M_PI*p.y)/std::sinh(M_PI);
    };

    std::vector<bool> is_bnd; build_boundary_mask(m, is_bnd);
    CSRPattern P; build_P1_CSRPattern(m, P);
    TArray<double> u(N), du(N,0.0), q(m.triangle_count()), rhs(N,0.0);
    TArray<double> u_tmp(N), r_cg(N), p_cg(N), Ap_cg(N);
    CSRMatrix S; init_csr_from_pattern(P, S);
    make_standard_init(m, f, u);

    MinimalGraphSolver helper(m, f);
    double area = helper.compute_denominator(q, u);

    constexpr double HUGE_DIAG=1e30, c_armijo=1e-4, rho=0.5;
    const double tolCG=1e-14;

    std::vector<std::string> header={"iter","kappa","residual_norm","correction_norm","area"};
    std::vector<std::vector<double>> rows;
    printf("%-6s %-14s %-14s %-14s %-14s\n","iter","kappa","residual","correction","area");

    for (size_t k=0; k<50; ++k) {
        build_P1_stiffness_matrix_NS(m,P,S,q.data,u.data,area);
        build_P1_rhs_NS(m,q.data,u.data,rhs);

        double res_sq=0;
        for (size_t i=0;i<N;++i) if (!is_bnd[i]) res_sq+=rhs[i]*rhs[i];
        double res_norm=std::sqrt(res_sq);

        /* estimate kappa before boundary enforcement */
        double kappa = estimate_condition_number(S, m, 50);

        for (size_t i=0;i<m.boundary.size;++i) {
            uint32_t bi=m.boundary[i]; S(bi,bi)=HUGE_DIAG; rhs[bi]=0.0;
        }
        memset(du.data,0,N*sizeof(double));
        double cg_err=0;
        conjugate_gradient_solve(S,rhs.data,du.data,r_cg.data,p_cg.data,Ap_cg.data,
                                 &cg_err,tolCG,10000,false);
        double corr_sq=blas_dot(du.data,du.data,N);

        double alpha=1.0; bool flag=true;
        while (flag) {
            for (size_t i=0;i<N;++i) u_tmp[i]=u[i]+alpha*du[i];
            double e_tmp=helper.compute_denominator(q,u_tmp);
            if (e_tmp<area-c_armijo*alpha*corr_sq||alpha<=1.0) flag=false;
            else alpha*=rho;
        }
        for (size_t i=0;i<N;++i) u[i]+=alpha*du[i];
        area=helper.compute_denominator(q,u);
        double corr_norm=std::sqrt(corr_sq);

        printf("%-6zu %-14.4e %-14.4e %-14.4e %-14.6f\n",k,kappa,res_norm,corr_norm,area);
        rows.push_back({(double)k,kappa,res_norm,corr_norm,area});
        if (corr_norm<1e-12) break;
    }
    write_csv("experiment_4_cond_newton_path.csv", header, rows);
}

/* =========================================================================
 * EXPERIMENT 5 — Condition number scaling with mesh refinement
 * ========================================================================= */
static void experiment_5()
{
    printf("\n=== Experiment 5: Condition Number Scaling with Mesh Refinement ===\n");

    auto f = [](const Vec2d &p) {
        return std::sin(M_PI*p.x)*std::sinh(M_PI*p.y)/std::sinh(M_PI);
    };

    std::vector<std::string> header={"level","N_subdiv","h","N_vtx","kappa"};
    std::vector<std::vector<double>> rows;
    printf("%-6s %-8s %-10s %-10s %-14s\n","level","N","h","N_vtx","kappa");

    for (size_t li=0; li<5; ++li) {
        size_t Nsub = 2u<<li;
        Mesh m; build_square_mesh(&m,Nsub,1.0,1.0);
        double h=1.0/(double)Nsub;

        MinimalGraphSolver solver(m,f);
        solver.clear_solution(true);
        solver.do_iterate_Newton(500,1e-12,1.0);

        size_t N=m.vertex_count();
        CSRPattern P; build_P1_CSRPattern(m,P);
        TArray<double> q(m.triangle_count());
        CSRMatrix S; init_csr_from_pattern(P,S);
        double area=solver.compute_denominator(q,solver.u);
        build_P1_stiffness_matrix_NS(m,P,S,q.data,solver.u.data,area);
        /* boundary enforcement before kappa (the new estimator projects it out) */
        for (size_t i=0;i<m.boundary.size;++i)
            S(m.boundary[i],m.boundary[i])=1e30;

        double kappa=estimate_condition_number(S,m,60);
        printf("%-6zu %-8zu %-10.4f %-10zu %-14.4e\n",li,Nsub,h,N,kappa);
        rows.push_back({(double)li,(double)Nsub,h,(double)N,kappa});
    }
    write_csv("experiment_5_cond_scaling.csv", header, rows);
}

/* =========================================================================
 * EXPERIMENT 6 — Effect of boundary data magnitude
 * ========================================================================= */
static void experiment_6()
{
    printf("\n=== Experiment 6: Effect of Boundary Data Magnitude ===\n");

    Mesh m; build_square_mesh(&m,10,1.0,1.0);
    size_t N=m.vertex_count();

    auto g_base=[](const Vec2d &p){
        return std::sin(M_PI*p.x)*std::sinh(M_PI*p.y)/std::sinh(M_PI);
    };

    std::vector<std::string> header={"lambda","converged","n_iter","final_area","kappa"};
    std::vector<std::vector<double>> rows;
    printf("%-8s %-12s %-8s %-14s %-14s\n","lambda","converged","n_iter","area","kappa");

    for (double lam : {0.5,1.0,2.0,4.0,8.0,16.0}) {
        auto g=[&](const Vec2d &p){ return lam*g_base(p); };
        TArray<double> u0(N); make_standard_init(m,g,u0);
        std::vector<NewtonRecord> recs;
        bool conv=run_newton_instrumented(m,g,u0,500,1e-12,recs,1.0);
        double fa=recs.empty()?0.0:recs.back().area;

        CSRPattern P; build_P1_CSRPattern(m,P);
        TArray<double> q(m.triangle_count());
        CSRMatrix S; init_csr_from_pattern(P,S);
        MinimalGraphSolver helper(m,g);
        double area=helper.compute_denominator(q,u0);
        build_P1_stiffness_matrix_NS(m,P,S,q.data,u0.data,area);
        for (size_t i=0;i<m.boundary.size;++i)
            S(m.boundary[i],m.boundary[i])=1e30;
        double kappa=estimate_condition_number(S,m,50);

        printf("lam=%5.1f  conv=%s  iter=%3zu  area=%.6f  kappa=%.4e\n",
               lam,conv?"yes":"no",recs.size(),fa,kappa);
        rows.push_back({lam,conv?1.0:0.0,(double)recs.size(),fa,kappa});
    }
    write_csv("experiment_6_boundary_magnitude.csv", header, rows);
}

/* =========================================================================
 * EXPERIMENT 7 — Energy decrease monitoring
 * ========================================================================= */
static void experiment_7()
{
    printf("\n=== Experiment 7: Energy Decrease Monitoring ===\n");

    Mesh m; build_square_mesh(&m,10,1.0,1.0);
    auto f=[](const Vec2d &p){
        return std::sin(M_PI*p.x)*std::sinh(M_PI*p.y)/std::sinh(M_PI);
    };
    TArray<double> u0(m.vertex_count()); make_standard_init(m,f,u0);

    std::vector<NewtonRecord> recs;
    bool conv=run_newton_instrumented(m,f,u0,300,1e-12,recs,1.0);
    printf("Converged: %s after %zu iterations\n",conv?"yes":"no",recs.size());

    bool monotone=true;
    for (size_t k=1;k<recs.size();++k)
        if (recs[k].area>recs[k-1].area+1e-12) monotone=false;
    printf("Energy monotonically decreasing: %s\n",monotone?"YES":"NO");

    std::vector<std::string> header={"iter","area","correction_norm","alpha","delta_area"};
    std::vector<std::vector<double>> rows;
    for (size_t k=0;k<recs.size();++k) {
        double da=(k>0)?recs[k].area-recs[k-1].area:0.0;
        printf("iter %3zu  area=%.8f  |du|=%.3e  alpha=%.4f  dA=%.3e\n",
               k,recs[k].area,recs[k].correction_norm,recs[k].alpha,da);
        rows.push_back({(double)k,recs[k].area,recs[k].correction_norm,recs[k].alpha,da});
    }
    write_csv("experiment_7_energy.csv", header, rows);
}

/* =========================================================================
 * EXPERIMENT 8 — Two sources of ill-conditioning (Scherk family)
 *
 * Test family: u_alpha(x,y) = (1/alpha)*ln(cos(alpha*x)/cos(alpha*y))
 * on [-1,1]^2, valid MSE solution for alpha in (0, pi/2).
 * M = tan(alpha).  Non-uniform gradient ensures d_K varies across elements.
 *
 * Part A: fix alpha=pi/4 (M=1), vary N  ->  kappa ~ h^{-2}
 * Part B: fix N=32, vary alpha           ->  kappa ~ (1+M^2)^{3/2}
 * Part C: C_hat = kappa*h^2/(1+M^2)^{3/2} should be ~constant
 * ========================================================================= */
static void experiment_8()
{
    printf("\n=== Experiment 8: Sources of Ill-Conditioning (Scherk) ===\n");
    printf("Theory: kappa(J) ~ C * h^{-2} * (1+M^2)^{3/2}\n\n");

    std::vector<std::string> hdr={
        "part","alpha","N","h","M",
        "kappa","log_slope","sqrt_kappa","cg_iters","cg_over_sqrtk"
    };
    std::vector<std::vector<double>> rows;

    /* Worker: converge Newton, build J at U*, estimate kappa and CG iters */
    auto run_one = [&](const Mesh &m, double alpha,
                       std::mt19937 &rng) -> std::pair<double,int>
    {
        const size_t N=m.vertex_count();
        auto u_a=[alpha](const Vec2d &p){
            return (1.0/alpha)*std::log(std::cos(alpha*p.x)/std::cos(alpha*p.y));
        };

        TArray<double> u(N); make_standard_init(m,u_a,u);
        { std::vector<NewtonRecord> recs;
          run_newton_instrumented(m,u_a,u,500,1e-12,recs,1.0); }

        CSRPattern P; build_P1_CSRPattern(m,P);
        TArray<double> q(m.triangle_count());
        MinimalGraphSolver helper(m,u_a);
        double area=helper.compute_denominator(q,u);
        CSRMatrix S; init_csr_from_pattern(P,S);
        build_P1_stiffness_matrix_NS(m,P,S,q.data,u.data,area);
        for (size_t i=0;i<m.boundary.size;++i)
            S(m.boundary[i],m.boundary[i])=1e30;

        double kappa=estimate_condition_number(S,m,80);

        std::vector<bool> is_bnd; build_boundary_mask(m,is_bnd);
        std::uniform_real_distribution<double> dist(-1.0,1.0);
        std::vector<double> b(N,0.0); double bnorm=0;
        for (size_t i=0;i<N;++i) if (!is_bnd[i]){ b[i]=dist(rng); bnorm+=b[i]*b[i]; }
        bnorm=std::sqrt(bnorm);
        for (size_t i=0;i<N;++i) b[i]/=bnorm;
        int cg=count_cg_iters(S,b,1e-6);
        return {kappa,cg};
    };

    std::mt19937 rng(0xC0FFEE42u);

    /* Part A */
    const double alpha_A=M_PI/4.0, M_A=std::tan(alpha_A);
    printf("Part A: Spatial penalty [alpha=pi/4, M=%.3f, vary N]\n",M_A);
    printf("%-5s %-8s %-12s %-7s %-10s %-7s %-8s\n",
           "N","h","kappa","slope","sqrt(k)","CG","CG/sqrtK");

    double prev_kA=0,prev_hA=0;
    for (int Ns : {4,8,16,32,64}) {
        Mesh m; build_square_mesh(&m,(size_t)Ns,1.0,1.0);
        rescale_and_recenter_mesh(m);
        double h=2.0/Ns;
        auto [kappa,cg]=run_one(m,alpha_A,rng);
        double slope=(prev_kA>0)?std::log(kappa/prev_kA)/std::log(h/prev_hA):0.0;
        double sqrtk=std::sqrt(kappa), ratio=cg/sqrtk;
        printf("%-5d %-8.4f %-12.3e %-7.2f %-10.1f %-7d %-8.3f\n",
               Ns,h,kappa,slope,sqrtk,cg,ratio);
        rows.push_back({1,alpha_A,(double)Ns,h,M_A,kappa,slope,sqrtk,(double)cg,ratio});
        prev_kA=kappa; prev_hA=h;
    }

    /* Part B */
    printf("\nPart B: Nonlinear penalty [N=32, vary alpha -> M]\n");
    printf("%-8s %-8s %-12s %-7s %-10s %-7s %-8s\n",
           "alpha/pi","M","kappa","slope","sqrt(k)","CG","CG/sqrtK");

    double prev_kB=0,prev_MB=0;
    for (double frac : {0.25,0.5,0.75,0.9,0.95}) {
        double alpha=frac*M_PI/2.0, M=std::tan(alpha);
        Mesh m; build_square_mesh(&m,32,1.0,1.0);
        rescale_and_recenter_mesh(m);
        double h=2.0/32;
        auto [kappa,cg]=run_one(m,alpha,rng);
        double slope=(prev_kB>0&&prev_MB>0)?
            std::log(kappa/prev_kB)/std::log((1+M*M)/(1+prev_MB*prev_MB)):0.0;
        double sqrtk=std::sqrt(kappa), ratio=cg/sqrtk;
        printf("%-8.3f %-8.4f %-12.3e %-7.2f %-10.1f %-7d %-8.3f\n",
               frac,M,kappa,slope,sqrtk,cg,ratio);
        rows.push_back({2,alpha,32,h,M,kappa,slope,sqrtk,(double)cg,ratio});
        prev_kB=kappa; prev_MB=M;
    }

    /* Part C */
    printf("\nPart C: C_hat = kappa * h^2 / (1+M^2)^{3/2}\n");
    printf("         ");
    for (int Ns : {8,16,32}) printf("  N=%-5d",Ns);
    printf("\n");

    for (double frac : {0.25,0.5,0.75,0.9}) {
        double alpha=frac*M_PI/2.0, M=std::tan(alpha);
        printf("a=%.2f*pi/2 (M=%.2f)  ",frac,M);
        for (int Ns : {8,16,32}) {
            Mesh m; build_square_mesh(&m,(size_t)Ns,1.0,1.0);
            rescale_and_recenter_mesh(m);
            double h=2.0/Ns;
            auto [kappa,cg]=run_one(m,alpha,rng);
            double chat=kappa*h*h/std::pow(1.0+M*M,1.5);
            printf("  %-9.4f",chat);
            rows.push_back({3,alpha,(double)Ns,h,M,kappa,chat,std::sqrt(kappa),(double)cg,0.0});
        }
        printf("\n");
    }

    write_csv("experiment_8_illconditioning.csv", hdr, rows);
}

/* =========================================================================
 * EXPERIMENT 9 — Preconditioner comparison (Jacobi, SSOR, IC(0))
 *
 * Setup: Scherk family at alpha=pi/4 (M=1) on [-1,1]^2, same as
 *        Experiment 8 Part A.  Newton converges to U*, then we solve
 *        the Jacobian system J du = F with four solvers:
 *
 *          CG        unpreconditioned              (baseline)
 *          PCG-J     Jacobi M = diag(J)
 *          PCG-SSOR  SSOR   M = (D+L)D^{-1}(D+L^T)
 *          PCG-IC    IC(0)  M = R^T R
 *
 * Reported per mesh level:
 *   h, N_vtx, kappa(J), kappa(M^{-1}J), iters_{CG,J,SSOR,IC},
 *   speedup_{J,SSOR,IC}  (= CG iters / PCG iters)
 *
 * Theory check:
 *   CG    iters ~ sqrt(kappa(J))         ~ h^{-1}
 *   PCG-J iters ~ sqrt(kappa(M_J^{-1}J)) < CG but h^{-1} scaling persists
 *   PCG-IC iters ~ sqrt(kappa(M_IC^{-1}J)) ~ h^{-1/2} (hopefully)
 * ========================================================================= */
static void experiment_9()
{
    printf("\n=== Experiment 9: Preconditioner Comparison ===\n");
    printf("Theory: CG iters ~ sqrt(kappa) ~ h^{-1};  "
           "IC should give ~ h^{-1/2}\n\n");

    /* ------------------------------------------------------------------ */
    /* PCG iteration counter — mirrors count_cg_iters but calls pcg_solve  */
    /* ------------------------------------------------------------------ */
    auto count_pcg_iters = [](const CSRMatrix &A,
                               const std::vector<double> &b,
                               auto &precond,
                               double tol = 1e-6,
                               int max_iter = 200000) -> int
    {
        size_t n = A.rows;
        std::vector<double> x(n, 0.0), r(b), z(n, 0.0), p(n, 0.0), Ap(n, 0.0);

        /* z_0 = M^{-1} r_0 */
        precond.apply(r.data(), z.data());
        p = z;

        double rz  = 0.0;
        for (size_t i = 0; i < n; ++i) rz += r[i] * z[i];
        double b2  = 0.0;
        for (size_t i = 0; i < n; ++i) b2 += b[i] * b[i];
        if (b2 < 1e-100) return 0;
        const double tol2 = tol * tol * b2;
        double r2  = 0.0;
        for (size_t i = 0; i < n; ++i) r2 += r[i] * r[i];

        for (int k = 0; k < max_iter; ++k) {
            if (r2 < tol2) return k;
            A.mvp(p.data(), Ap.data());
            double pAp = 0.0;
            for (size_t i = 0; i < n; ++i) pAp += p[i] * Ap[i];
            if (std::fabs(pAp) < 1e-100) return k;
            double alpha = rz / pAp;
            for (size_t i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }
            precond.apply(r.data(), z.data());
            double rz_new = 0.0;
            for (size_t i = 0; i < n; ++i) rz_new += r[i] * z[i];
            double beta = rz_new / rz;
            for (size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
            r2  = 0.0;
            for (size_t i = 0; i < n; ++i) r2 += r[i] * r[i];
            rz  = rz_new;
        }
        return max_iter;
    };

    /* ------------------------------------------------------------------ */
    /* Condition number of M^{-1}A via power + inverse-power iteration     */
    /* ------------------------------------------------------------------ */
    auto estimate_precond_kappa = [](const CSRMatrix &A,
                                     auto &precond,
                                     const Mesh &m,
                                     size_t N_pow = 60) -> double
    {
        size_t N = m.vertex_count();
        std::vector<bool> is_bnd;
        build_boundary_mask(m, is_bnd);

        std::mt19937 rng(54321);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        TArray<double> v(N, 0.0), w(N, 0.0), z(N, 0.0);

        /* lambda_max of M^{-1}A via power iteration on v -> M^{-1}(A v) */
        for (size_t i = 0; i < N; ++i) v[i] = is_bnd[i] ? 0.0 : dist(rng);
        double vnm = std::sqrt(blas_dot(v.data, v.data, N));
        if (vnm > 0) for (size_t i = 0; i < N; ++i) v[i] /= vnm;

        double lambda_max = 1.0;
        for (size_t k = 0; k < N_pow; ++k) {
            A.mvp(v.data, w.data);
            for (size_t i = 0; i < N; ++i) if (is_bnd[i]) w[i] = 0.0;
            precond.apply(w.data, z.data);        /* z = M^{-1} (A v) */
            for (size_t i = 0; i < N; ++i) if (is_bnd[i]) z[i] = 0.0;
            lambda_max = blas_dot(v.data, z.data, N);
            double zn = std::sqrt(blas_dot(z.data, z.data, N));
            if (zn < 1e-30) break;
            for (size_t i = 0; i < N; ++i) v[i] = z[i] / zn;
        }

        /* lambda_min via inverse power: solve M^{-1}A y = v repeatedly */
        for (size_t i = 0; i < N; ++i) v[i] = is_bnd[i] ? 0.0 : dist(rng);
        vnm = std::sqrt(blas_dot(v.data, v.data, N));
        if (vnm > 0) for (size_t i = 0; i < N; ++i) v[i] /= vnm;

        TArray<double> y(N, 0.0), r_cg(N), p_cg(N), Ap_cg(N);
        double lambda_min = lambda_max;

        for (size_t inv_k = 0; inv_k < 5; ++inv_k) {
            /* Solve M^{-1}A y = v  by applying (M^{-1}A)^{-1} = A^{-1}M:
             * first apply M (already have v), then solve A y = M v via CG. */
            std::vector<double> Mv(N);
            /* M v: since M is SPD, apply M by "un-applying" M^{-1}:
             * we can't apply M directly from the precond interface.
             * Instead use Rayleigh quotient on (A v):  lambda_min ~ v^T A v */
            A.mvp(v.data, w.data);
            for (size_t i = 0; i < N; ++i) if (is_bnd[i]) w[i] = 0.0;
            lambda_min = blas_dot(v.data, w.data, N);

            /* solve A y = v */
            memset(y.data, 0, N * sizeof(double));
            double cg_err = 0.0;
            conjugate_gradient_solve(A, v.data, y.data,
                                     r_cg.data, p_cg.data, Ap_cg.data,
                                     &cg_err, 1e-10, 5000, false);
            for (size_t i = 0; i < N; ++i) if (is_bnd[i]) y[i] = 0.0;

            /* apply M^{-1} to get next iterate of inverse power */
            precond.apply(y.data, z.data);
            for (size_t i = 0; i < N; ++i) if (is_bnd[i]) z[i] = 0.0;

            double zn = std::sqrt(blas_dot(z.data, z.data, N));
            if (zn < 1e-30) break;
            for (size_t i = 0; i < N; ++i) v[i] = z[i] / zn;
        }

        if (lambda_min < 1e-30) return 1e30;
        return lambda_max / lambda_min;
    };

    /* ------------------------------------------------------------------ */
    /* Main loop over mesh levels                                           */
    /* ------------------------------------------------------------------ */
    std::vector<std::string> hdr = {
        "N", "h", "N_vtx",
        "kappa_J",
        "kappa_precond_jacobi", "kappa_precond_ssor", "kappa_precond_ic",
        "iters_cg", "iters_jacobi", "iters_ssor", "iters_ic",
        "speedup_jacobi", "speedup_ssor", "speedup_ic"
    };
    std::vector<std::vector<double>> rows;

    printf("%-5s %-8s %-8s %-12s  %-12s %-12s %-12s  "
           "%-7s %-7s %-7s %-7s  %-8s %-8s %-8s\n",
           "N", "h", "N_vtx", "kappa(J)",
           "k(M_J^-1J)", "k(M_S^-1J)", "k(M_IC^-1J)",
           "CG", "J", "SSOR", "IC",
           "spdup_J", "spdup_S", "spdup_IC");

    const double alpha = M_PI / 4.0;
    auto u_scherk_alpha = [alpha](const Vec2d &p) {
        return (1.0 / alpha) *
               std::log(std::cos(alpha * p.x) / std::cos(alpha * p.y));
    };

    std::mt19937 rng(0xC0FFEE42u);

    for (int Ns : {4, 8, 16, 32, 64}) {
        Mesh m;
        build_square_mesh(&m, (size_t)Ns, 1.0, 1.0);
        rescale_and_recenter_mesh(m);
        const double h    = 2.0 / Ns;
        const size_t Nvtx = m.vertex_count();

        /* -- Converge Newton to U* -- */
        TArray<double> u(Nvtx);
        make_standard_init(m, u_scherk_alpha, u);
        {
            std::vector<NewtonRecord> recs;
            run_newton_instrumented(m, u_scherk_alpha, u, 500, 1e-12, recs, 1.0);
        }

        /* -- Build Jacobian at U* -- */
        CSRPattern P;
        build_P1_CSRPattern(m, P);
        TArray<double> q(m.triangle_count());
        MinimalGraphSolver helper(m, u_scherk_alpha);
        double area = helper.compute_denominator(q, u);
        CSRMatrix J;
        init_csr_from_pattern(P, J);
        build_P1_stiffness_matrix_NS(m, P, J, q.data, u.data, area);
        for (size_t i = 0; i < m.boundary.size; ++i)
            J(m.boundary[i], m.boundary[i]) = 1e30;

        /* -- Random interior RHS -- */
        std::vector<bool> is_bnd;
        build_boundary_mask(m, is_bnd);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<double> b(Nvtx, 0.0);
        double bnorm = 0.0;
        for (size_t i = 0; i < Nvtx; ++i)
            if (!is_bnd[i]) { b[i] = dist(rng); bnorm += b[i] * b[i]; }
        bnorm = std::sqrt(bnorm);
        for (size_t i = 0; i < Nvtx; ++i) b[i] /= bnorm;

        /* -- kappa(J) -- */
        double kappa_J = estimate_condition_number(J, m, 80);

        /* -- Build all three preconditioners -- */
        JacobiPrecond<CSRMatrix> pJ;   pJ.build(J);
        SSORPrecond<CSRMatrix>   pS;   pS.build(J);
        IncompleteCholeskyPrecond<CSRMatrix> pIC; pIC.build(J);

        /* -- kappa(M^{-1}J) for each preconditioner -- */
        double kappa_pJ  = estimate_precond_kappa(J, pJ,  m, 80);
        double kappa_pS  = estimate_precond_kappa(J, pS,  m, 80);
        double kappa_pIC = estimate_precond_kappa(J, pIC, m, 80);

        /* -- Iteration counts -- */
        int iters_cg   = count_cg_iters(J, b, 1e-6);
        int iters_J    = count_pcg_iters(J, b, pJ,  1e-6);
        int iters_S    = count_pcg_iters(J, b, pS,  1e-6);
        int iters_IC   = count_pcg_iters(J, b, pIC, 1e-6);

        double spdup_J  = (iters_J  > 0) ? (double)iters_cg / iters_J  : 0.0;
        double spdup_S  = (iters_S  > 0) ? (double)iters_cg / iters_S  : 0.0;
        double spdup_IC = (iters_IC > 0) ? (double)iters_cg / iters_IC : 0.0;

        printf("%-5d %-8.4f %-8zu %-12.3e  %-12.3e %-12.3e %-12.3e  "
               "%-7d %-7d %-7d %-7d  %-8.2f %-8.2f %-8.2f\n",
               Ns, h, Nvtx, kappa_J,
               kappa_pJ, kappa_pS, kappa_pIC,
               iters_cg, iters_J, iters_S, iters_IC,
               spdup_J, spdup_S, spdup_IC);

        rows.push_back({
            (double)Ns, h, (double)Nvtx,
            kappa_J,
            kappa_pJ, kappa_pS, kappa_pIC,
            (double)iters_cg, (double)iters_J, (double)iters_S, (double)iters_IC,
            spdup_J, spdup_S, spdup_IC
        });

        pJ.free_data();
        pS.free_data();
        pIC.free_data();
    }

    write_csv("experiment_9_preconditioners.csv", hdr, rows);
}

/* =========================================================================
 * main
 * ========================================================================= */
int main(int argc, char **argv)
{
    int experiment=0;
    for (int i=1;i<argc-1;++i)
        if (std::string(argv[i])=="--experiment")
            experiment=std::atoi(argv[i+1]);

    if (experiment==0) {
        printf("Usage: %s --experiment N  (N = 1..8)\n", argv[0]);
        return EXIT_FAILURE;
    }
    switch (experiment) {
    case 1: experiment_1(); break;
    case 2: experiment_2(); break;
    case 3: experiment_3(); break;
    case 4: experiment_4(); break;
    case 5: experiment_5(); break;
    case 6: experiment_6(); break;
    case 7: experiment_7(); break;
    case 8: experiment_8(); break;
    case 9: experiment_9(); break;
    default:
        fprintf(stderr,"Unknown experiment %d (must be 1–9)\n",experiment);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}