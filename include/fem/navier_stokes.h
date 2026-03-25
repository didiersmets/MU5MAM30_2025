#include <vector>

#include "array.h"
#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif
#include "mesh.h"
#include "cholesky.h"

enum class SolverType { CG, CHOLESKY };

struct NavierStokesSolver {
	NavierStokesSolver(const Mesh &m, int degre = 1, SolverType solver = SolverType::CG, double nu = 1e-3, double dt = 0.002);/*ajout du degre pour P2*/
	const Mesh &m;
	int degre;
	size_t N;   // DoF
	std::vector<uint32_t> edge_ddls;
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

	// Cholesky solver fields
	SolverType solver_type = SolverType::CG;
	bool cholesky_ready = false;
	CSRPattern cholesky_A_pattern;  // Pattern for matrix A = M + nu*dt*S
	TArray<uint32_t> cholesky_perm;
	TArray<uint32_t> cholesky_iperm;
	CSRPattern cholesky_L_pattern;
	CSRMatrix cholesky_L;

	bool inited; // Initialization computes first residue and error

	size_t iter_max = 500;
	double tol = 1e-6;

	double t;

	void set_zero_mean(double *V);
	size_t compute_stream_function();
	void compute_transport(double *T);
	void time_step(double dt, double nu);
};
