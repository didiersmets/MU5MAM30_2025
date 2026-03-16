#include "array.h"
#include <vector>
#include <stdint.h>

#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
    #include "fem_matrix.h"
#else
    #include "sparse_matrix.h"
    #include "sparse_cholesky.h"  
#endif
#include "mesh.h"

struct NavierStokesSolver {
    NavierStokesSolver(const Mesh &m);
    const Mesh &m;
    size_t N;   // DoF
    double vol; // Surface(m), used for insuring zero mean to omega and psi

    TArray<double> omega;
    TArray<double> Momega;
    TArray<double> psi;

#if USE_FEM_MATRIX
    FEMatrix S; // Stiffness matrix
    FEMatrix M; // Mass matrix
#else
    CSRPattern P; // Pattern arrays
    CSRMatrix S;  // Stiffness matrix
    CSRMatrix M;  // Mass matrix

    // Solver Cholesky
    bool chol_ready = false;
    double chol_alpha = -1.0;              // alpha = nu * dt used to build A

    CSRMatrix A;                           
    std::vector<uint32_t> etree_parent;    // elimination tree parent
    CSRPattern PL;                         // pattern of L
    CSRMatrix L;                           // Cholesky factor (shares pattern with PL)
    std::vector<std::vector<ColAdjEntry>> col_adj; // for L^T solve
    TArray<double> ysolve;                 // temporary buffer y (Ly=b)
#endif

    TArray<double> r;  // current residue r = Mf - Su
    TArray<double> p;  // internal for cg (reused as RHS buffer for omega solve)
    TArray<double> Ap; // internal for cg (unused for omega once direct solver is used)

    bool inited; // Initialization computes first residue and error

    size_t iter_max = 500;
    double tol = 1e-6;

    double t;

    void set_zero_mean(double *V);
    size_t compute_stream_function();
    void compute_transport(double *T);
    void time_step(double dt, double nu);

#if !USE_FEM_MATRIX
    void init_direct_solver(double dt, double nu);
#endif
};
