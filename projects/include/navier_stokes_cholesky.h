#include "array.h"
#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif
#include "mesh.h"

/* This is basically identical to navier_stokes provided in class, with minor changes to adapt
to Cholesky setting. */

struct NavierStokesSolverCholesky {
	NavierStokesSolverCholesky(const Mesh &m, double dt, double nu);
	const Mesh &m;
	size_t N;   // DoF
	double vol; // Surface(m), used for insuring zero mean to omega and psi

	TArray<double> omega;
	TArray<double> Momega;
	TArray<double> psi;
	TArray<double> transport;
#if USE_FEM_MATRIX
	FEMatrix S; // Stiffness matrix
	FEMatrix M; // Mass matrix
#else
	CSRPattern P; // Pattern arrays
	CSRMatrix S;  // Stiffness matrix
	CSRMatrix M;  // Mass matrix
#endif
	TArray<uint32_t> parent;  
	CSRPattern L_pattern ;  
	CSRMatrix L_S;
	CSRMatrix L_SM;

	double t;
	double dt;
	double nu;

	void set_zero_mean(double *V);
	void compute_stream_function();
	void compute_transport(double *T);
	void time_step();
};
