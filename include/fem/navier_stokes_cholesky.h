#include "array.h"
#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif

#include "cholesky.h"   /* CholeskySolver */
#include "mesh.h"
#include <string.h>

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
#endif

	TArray<double> r;  // scratch
	TArray<double> p;  // scratch
	TArray<double> Ap; // scratch

	/* Cholesky factorization of S (for stream function solve).
	 * NOTE: S has a 1D kernel (constants); we handle this by
	 * using set_zero_mean before and after the solve, exactly
	 * as we did with CG.  The factorization itself ignores the
	 * kernel — it will produce a valid solution modulo a constant,
	 * which set_zero_mean then removes.                              */
	CholeskySolver chol_S;

	/* Cholesky factorization of  A_dt = M + nu*dt*S.
	 * Must be rebuilt whenever dt or nu changes via setup_cholesky(). */
	CholeskySolver chol_A;

	/* CSRMatrix wrapper for  M + nu*dt*S  (owns its data array).
	 * row_start and col are shared with P (same sparsity pattern).  */
	CSRMatrix A_dt;   /* M + nu*dt*S */
	double    last_dt = -1.0;
	double    last_nu = -1.0;

	bool inited;

	size_t iter_max = 500;
	double tol = 1e-6;

	double t;

	/* Call this once (or whenever dt/nu change) before starting the loop. */
	void setup_cholesky(double dt, double nu);

	void set_zero_mean(double *V);
	size_t compute_stream_function();
	void compute_transport(double *T);
	void time_step(double dt, double nu);
};