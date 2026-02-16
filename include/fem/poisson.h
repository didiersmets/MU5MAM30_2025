#pragma once

#include "array.h"
#include "mesh.h"
#include "fem_type.h"   // <--- NUEVO

#define USE_FEM_MATRIX false // Lo cambie
#if USE_FEM_MATRIX
#include "fem_matrix.h"
#else
#include "sparse_matrix.h"
#endif

struct PoissonSolver {
    // -------------------------
    // NUEVO CONSTRUCTOR
    // -------------------------
    PoissonSolver(const Mesh &m, FEMType fem_type);

    const Mesh &m;
    FEMType fem_type;   // <--- NUEVO

    size_t N; // DoF
    double vol; // Surface(m), used for insuring zero mean to f and u

    TArray<double> f;
    TArray<double> u;

    // --- Para condiciones de contorno Dirichlet --- (para estudiar el error)
    TArray<char> is_dirichlet;     // 1 si es nodo Dirichlet
    TArray<double> u_dirichlet;    // valor Dirichlet
    void apply_dirichlet();        // función que modifica A y f


#if USE_FEM_MATRIX
    FEMatrix A; // Stiffness matrix
    FEMatrix M; // Mass matrix
#else
    CSRPattern P; // Pattern arrays
    CSRMatrix A; // Stiffness matrix
    CSRMatrix M; // Mass matrix
#endif

    TArray<double> r;  // current residue r = Mf - Su
    TArray<double> p;  // internal for cg
    TArray<double> Ap; // internal for cg

    bool inited;
    size_t iterate;
    double b2;
    double r2;
    bool converged;
    double rel_error;

    void clear_solution();
    void init_cg();
    void set_zero_mean(double *V);
    void do_iterate(size_t max_iter, double tol);
};
