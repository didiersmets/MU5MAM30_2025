#include <assert.h>
#include <cmath>
#include <cstring>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

struct V3 {
    double x, y, z;
};

static inline V3 v3(double x, double y, double z) { return {x,y,z}; }

static inline V3 add(const V3& a, const V3& b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline V3 sub(const V3& a, const V3& b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline V3 mul(double s, const V3& a)     { return {s*a.x, s*a.y, s*a.z}; }

static inline double dot(const V3& a, const V3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

static inline V3 cross(const V3& a, const V3& b) {
    return { a.y*b.z - a.z*b.y,
             a.z*b.x - a.x*b.z,
             a.x*b.y - a.y*b.x };
}

// Gradients of barycentric basis on a 3D-embedded triangle:
// n = (B-A) x (C-A), nn = |n|^2
// grad(lambda_A) = (n x (B-C)) / nn, ...
static inline void tri_grad_phi(const V3& A, const V3& B, const V3& C,
                                V3& gA, V3& gB, V3& gC,
                                V3& n, double& area)
{
    const V3 AB = sub(B, A);
    const V3 AC = sub(C, A);
    n = cross(AB, AC);
    const double nn = dot(n, n);
    if (nn <= 0.0) { area = 0.0; gA=gB=gC=v3(0,0,0); return; }

    area = 0.5 * std::sqrt(nn);

    // gA = (n x (B-C)) / nn, gB=..., gC=...
    gA = mul(1.0/nn, cross(n, sub(B, C)));
    gB = mul(1.0/nn, cross(n, sub(C, A)));
    gC = mul(1.0/nn, cross(n, sub(A, B)));
}

// ------------------------------------------------------------------

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
    // mean(V) := (1/vol) * 1^T M V
    // Implemented exactly like your Poisson: Ap = M V, s = sum(Ap), V -= s/vol
    M.mvp(V, Ap.data);
    const double s = blas_sum_in_place(Ap.data, N);
    const double c = s / vol;
    for (size_t i = 0; i < N; ++i) {
        V[i] -= c;
    }
}

size_t NavierStokesSolver::compute_stream_function()
{
    // Solve: S psi = M omega   
    // We'll do a simple CG loop.

    // rhs = M*omega into r
    M.mvp(omega.data, r.data);

    // Ensure compatibility (important on closed surfaces)
    set_zero_mean(omega.data);      // omega must have zero mean
    M.mvp(omega.data, r.data);      // refresh rhs after projection

    // Project initial guess psi to zero mean (helps CG stay in correct subspace)
    set_zero_mean(psi.data);

    // r = rhs - S*psi
    S.mvp(psi.data, Ap.data);
    blas_axpy(-1.0, Ap.data, r.data, N);

    const double b2 = blas_dot(r.data, r.data, N); // use rhs-norm proxy
    if (b2 == 0.0) return 0;

    blas_copy(r.data, p.data, N);

    double r2 = blas_dot(r.data, r.data, N);
    double rel = std::sqrt(r2 / b2);

    size_t iter = 0;
    while (iter < iter_max && rel > tol) {
        S.mvp(p.data, Ap.data);
        const double pAp = blas_dot(p.data, Ap.data, N);
        if (pAp <= 0.0) break; // should not happen in zero-mean subspace

        const double alpha = r2 / pAp;

        // psi += alpha * p
        blas_axpy(alpha, p.data, psi.data, N);

        // r -= alpha * Ap
        blas_axpy(-alpha, Ap.data, r.data, N);

        const double r2_new = blas_dot(r.data, r.data, N);
        const double beta = r2_new / r2;

        // p = r + beta * p
        for (size_t i = 0; i < N; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        r2 = r2_new;
        rel = std::sqrt(r2 / b2);
        ++iter;
    }

    // Clean any drift in nullspace
    set_zero_mean(psi.data);

    return iter;
}

void NavierStokesSolver::compute_transport(double *T)
{
    // Conservative-form transport:

    std::memset(T, 0, N * sizeof(double));

    const TArray<uint32_t> &idx = m.indices;
    const size_t nt = m.triangle_count();

    for (size_t t = 0; t < nt; ++t) {
        const uint32_t ia = idx[3*t + 0];
        const uint32_t ib = idx[3*t + 1];
        const uint32_t ic = idx[3*t + 2];

        const auto &Af = m.positions[ia];
        const auto &Bf = m.positions[ib];
        const auto &Cf = m.positions[ic];

        const V3 A = v3(Af.x, Af.y, Af.z);
        const V3 B = v3(Bf.x, Bf.y, Bf.z);
        const V3 C = v3(Cf.x, Cf.y, Cf.z);

        V3 gA, gB, gC, n;
        double area = 0.0;
        tri_grad_phi(A, B, C, gA, gB, gC, n, area);
        if (area == 0.0) continue;

        const double nn = dot(n, n);
        const double inv_norm_n = 1.0 / std::sqrt(nn);
        const V3 nhat = mul(inv_norm_n, n);

        // grad(psi) building 
        const V3 grad_psi = add(add(mul(psi[ia], gA), mul(psi[ib], gB)), mul(psi[ic], gC));

        // u = n_hat x ∇ψ
        const V3 u = cross(nhat, grad_psi);

        // ω average on triangle
        const double wavg = (omega[ia] + omega[ib] + omega[ic]) / 3.0;

        // T_i += area * wavg * (u · grad(psi_i))
        T[ia] += area * wavg * dot(u, gA);
        T[ib] += area * wavg * dot(u, gB);
        T[ic] += area * wavg * dot(u, gC);
    }

    // Optional: keep transport mean-free (helps keep omega mean-free without relying on post-projection)
    // set_zero_mean(T);
}

void NavierStokesSolver::time_step(double dt, double nu)
{
    compute_stream_function();

    /**********************************************************************
     * Solve:
     *  (M + nu*dt*S) omega^{n+1} = M*omega^n + dt*T(omega,psi)^n
     *********************************************************************/

    // T stored into r (reuse memory)
    compute_transport(r.data);

    // b = M*omega + dt*T   (store b in Momega)
    M.mvp(omega.data, Momega.data);
    blas_axpy(dt, r.data, Momega.data, N);

    // CG to solve A x = b with A = M + alpha*S, alpha = nu*dt
    const double alphaS = nu * dt;

    // initial guess x = omega (warm start)
    // residual r = b - A*x
    // Ap used as temp for A*x
    // A*x = M*x + alphaS * (S*x)
    M.mvp(omega.data, Ap.data);
    S.mvp(omega.data, p.data);                 // reuse p as temp: p = S*x
    blas_axpy(alphaS, p.data, Ap.data, N);     // Ap = Mx + alphaS*Sx

    // r = b - Ap
    blas_copy(Momega.data, r.data, N);
    blas_axpy(-1.0, Ap.data, r.data, N);

    const double b2 = blas_dot(Momega.data, Momega.data, N);
    if (b2 == 0.0) {
        set_zero_mean(omega.data);
        t += dt;
        return;
    }

    blas_copy(r.data, p.data, N);

    double r2 = blas_dot(r.data, r.data, N);
    double rel = std::sqrt(r2 / b2);

    size_t iter = 0;
    while (iter < iter_max && rel > tol) {
        // Ap = A*p = M*p + alphaS*S*p
        M.mvp(p.data, Ap.data);
        S.mvp(p.data, Momega.data);                 // reuse Momega as temp: Sp
        blas_axpy(alphaS, Momega.data, Ap.data, N); // Ap = Mp + alphaS*Sp

        const double pAp = blas_dot(p.data, Ap.data, N);
        if (pAp <= 0.0) break; // should be SPD because M is SPD and alphaS>=0

        const double alpha = r2 / pAp;

        // omega += alpha*p
        blas_axpy(alpha, p.data, omega.data, N);

        // r -= alpha*Ap
        blas_axpy(-alpha, Ap.data, r.data, N);

        const double r2_new = blas_dot(r.data, r.data, N);
        const double beta = r2_new / r2;

        // p = r + beta*p
        for (size_t i = 0; i < N; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        r2 = r2_new;
        rel = std::sqrt(r2 / b2);
        ++iter;
    }

    // mean-free vorticity
    set_zero_mean(omega.data);

    t += dt;
}