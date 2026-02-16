#include "array.h"
#define USE_FEM_MATRIX false

#if USE_FEM_MATRIX
    #include "fem_matrix.h"
#else
    #include "sparse_matrix.h"
#endif

#include "mesh.h"
#include "fem_matrix.h"  // ✅ para FEMType (como en Poisson)
#include "fem_type.h"

struct NavierStokesSolver {
    // ✅ ahora recibe también el tipo de FEM (P1 o P2)
    NavierStokesSolver(const Mesh &m, FEMType fem_type);

    const Mesh &m;
    size_t N;   // DoF totales (vértices [+ aristas en P2])
    double vol; // Surface(m), usada para imponer media cero a omega y psi

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
#endif

    TArray<double> r;  // current residue r = Mf - Su
    TArray<double> p;  // internal for cg
    TArray<double> Ap; // internal for cg

    bool inited; // Initialization computes first residue and error

    size_t iter_max = 500;
    double tol = 1e-6;

    double t;

    void set_zero_mean(double *V);
    size_t compute_stream_function();
    void compute_transport(double *T);
    void time_step(double dt, double nu);
};
