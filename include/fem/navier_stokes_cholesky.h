#include "array.h"
#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif

#include "cholesky.h"
#include "mesh.h"
#include <string.h>

struct NavierStokesSolver {
	NavierStokesSolver(const Mesh &m, double dt, double nu); // dt and nu now required at construction
	const Mesh &m;
	size_t N;
	double vol;

	TArray<double> omega;
	TArray<double> Momega;
	TArray<double> psi;
#if USE_FEM_MATRIX
	FEMatrix S;
	FEMatrix M;
#else
	CSRPattern P;
	CSRMatrix S;
	CSRMatrix M;
#endif

	TArray<double> r;
	TArray<double> p;
	TArray<double> Ap;

	/* Cholesky factorization of S (for stream function solve) */
	CholeskySolver chol_S;

	/* Cholesky factorization of A = M + nu*dt*S (for time step solve) */
	CholeskySolver chol_A;

	/* CSRMatrix for M + nu*dt*S */
	CSRMatrix A_dt;

	bool inited;

	size_t iter_max = 500;
	double tol = 1e-6;

	double t;

	void set_zero_mean(double *V);
	size_t compute_stream_function();
	void compute_transport(double *T);
	void time_step(double dt, double nu);
};