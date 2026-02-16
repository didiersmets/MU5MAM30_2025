#include "poisson.h"

#include "P1.h"
#include "P2.h"          // <--- NUEVO
#include "array.h"
#include "conjugate_gradient.h"
#if USE_FEM_MATRIX
    #include "fem_matrix.h"
#else
    #include "sparse_matrix.h"
#endif
#include "mesh.h"
#include "tiny_blas.h"

PoissonSolver::PoissonSolver(const Mesh &m, FEMType fem_type)
    : m(m),
      fem_type(fem_type),
      N(0),
      vol(0.0),
      f(),
      u(),
      r(),
      p(),
      Ap()
{
#if USE_FEM_MATRIX
    // Versión FEMatrix
    if (fem_type == FEMType::P1) {
        build_P1_mass_matrix(m, M);
        build_P1_stiffness_matrix(m, A);
    } else { // P2
        build_P2_mass_matrix(m, M);
        build_P2_stiffness_matrix(m, A);
    }

    N = M.rows;  // número de grados de libertad
#else
    // Versión CSR
    if (fem_type == FEMType::P1) {
        build_P1_CSRPattern(m, P);
        build_P1_mass_matrix(m, P, M);
        build_P1_stiffness_matrix(m, P, A);
    } else { // P2
        build_P2_CSRPattern(m, P);
        build_P2_mass_matrix(m, P, M);
        build_P2_stiffness_matrix(m, P, A);
    }

    N = M.rows;  // número de grados de libertad
#endif

    // Redimensionar vectores según N
    f.resize(N);
    u.resize(N);
    r.resize(N);
    p.resize(N);
    Ap.resize(N);
    is_dirichlet.resize(N);
    u_dirichlet.resize(N);
    for (size_t i = 0; i < N; ++i) {
        is_dirichlet[i] = 0;
        u_dirichlet[i] = 0.0;
    }


    // Inicializar solución a cero
    for (size_t i = 0; i < N; ++i)
        u[i] = 0.0;

    vol = M.sum();
    inited = false;
    iterate = 0;
    converged = false;
}

void PoissonSolver::clear_solution()
{
    for (size_t i = 0; i < N; i++)
        u[i] = 0.0;

    inited = false;
    iterate = 0;
    converged = false;
}

void PoissonSolver::apply_dirichlet()
{
    // Para cada dof marcado como Dirichlet
    for (size_t i = 0; i < N; ++i) {
        if (!is_dirichlet[i]) continue;

        // 1. Anular FILA i de A
        for (uint32_t k = A.row_start[i]; k < A.row_start[i+1]; ++k) {
            A.data[k] = 0.0;
        }

        // 2. Poner 1 en la diagonal A(i,i)
        for (uint32_t k = A.row_start[i]; k < A.row_start[i+1]; ++k) {
            if (A.col[k] == i) {
                A.data[k] = 1.0;
                break;
            }
        }

        // 3. Anular COLUMNA i de A en el resto de filas
        for (size_t row = 0; row < N; ++row) {
            if (row == i) continue;

            for (uint32_t k = A.row_start[row]; k < A.row_start[row+1]; ++k) {
                if (A.col[k] == i) {
                    A.data[k] = 0.0;
                }
            }
        }

        // 4. RHS = valor Dirichlet
        f[i] = u_dirichlet[i];
    }
}



void PoissonSolver::init_cg()
{
    apply_dirichlet();

    double *F  = f.data;
    double *U  = u.data;
    double *R  = r.data;
    double *P  = p.data;
    double *AP = Ap.data;

    // 1) Construir RHS
    if (fem_type == FEMType::P1) {
        // P1: f = f(x_i), necesitamos M*f
        M.mvp(F, R);          // R = b = M*f
    } else {
        // P2: f ya es b = ∫ f φ_i
        blas_copy(F, R, N);   // R = b
    }

    // 2) r0 = b - A*u0
    A.mvp(U, AP);
    blas_axpy(-1.0, AP, R, N);   // R = R - A*u0

    // 3) p0 = r0
    blas_copy(R, P, N);

    // 4) normas
    r2 = blas_dot(R, R, N);
    b2 = blas_dot(R, R, N);   // o similar
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
