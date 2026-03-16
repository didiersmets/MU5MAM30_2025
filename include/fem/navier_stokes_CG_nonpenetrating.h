#pragma once

#include "array.h"
#define USE_FEM_MATRIX false
#if USE_FEM_MATRIX
	#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif

#include "mesh.h"
#include <string.h>
#include <vector>

struct NavierStokesSolver {
	NavierStokesSolver(const Mesh &m);
	const Mesh &m;
	size_t N;
	double vol;

	std::vector<bool> is_boundary;
	bool has_boundary;

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

	bool inited;

	size_t iter_max = 500;
	double tol = 1e-6;

	double t;

	void set_zero_mean(double *V);
	size_t compute_stream_function();
	void compute_transport(double *T);
	void time_step(double dt, double nu);
};
