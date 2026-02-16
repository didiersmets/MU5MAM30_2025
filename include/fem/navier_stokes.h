#include "array.h"
#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
#include "fem_matrix.h"
#else
#include "sparse_matrix.h"
#endif
#include "mesh.h"

struct NavierStokesSolver
{
	NavierStokesSolver(Mesh &m, const bool &use_fem_P2);
	Mesh &m;
	size_t N;	// DoF
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
#endif
	TArray<double> r;  // current residue r = Mf - Su
	TArray<double> p;  // internal for cg
	TArray<double> Ap; // internal for cg

	bool inited; // Initialization computes first residue and error

	size_t iter_max = 500;
	double tol = 1e-6;
	double *rel_error; // sqrt(r2 / b2)
	double b2;		   // ||Mf||^2
	double r2;		   // current ||r||^2

	double t;

	bool use_fem_P2;

	void set_zero_mean(double *V);
	size_t compute_stream_function();
	void compute_transport(double *T);
	void time_step(double dt, double nu);
};
