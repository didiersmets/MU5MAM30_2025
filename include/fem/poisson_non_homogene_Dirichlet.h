#pragma once

#include "array.h"
#include <functional>
#include <vector>

#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
#include "fem_matrix.h"
#else
#include "sparse_matrix.h"
#endif
#include "mesh.h"

/* Type of the boundary value function:
 * takes (x, y, z) coordinates of a vertex, returns the Dirichlet value g. */
using BoundaryFunc = std::function<double(double x, double y, double z)>;

struct PoissonSolverDirichlet {
	PoissonSolverDirichlet(const Mesh &m, BoundaryFunc g);

	const Mesh &m;
	size_t N;    // DoF
	double vol;  // Surface area

	/* Boundary info */
	std::vector<bool> is_boundary;
	bool has_boundary;
	TArray<double> g_values; // g evaluated at each boundary vertex

	TArray<double> f;  // RHS (source term)
	TArray<double> u;  // solution
#if USE_FEM_MATRIX
	FEMatrix A;
	FEMatrix M;
#else
	CSRPattern P;
	CSRMatrix A;
	CSRMatrix M;
#endif
	TArray<double> r;
	TArray<double> p;
	TArray<double> Ap;

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

private:
	void apply_dirichlet_to_rhs();
};
