#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <functional>

#include "minimal_graph.h"

#include "P1.h"
#include "conjugate_gradient.h"
#include "tiny_blas.h"

namespace {

constexpr double HUGE_DIAG = 1e30;

inline Vec2d xy(const Mesh &m, uint32_t i)
{
	return Vec2d(static_cast<double>(m.positions[i].x),
		     static_cast<double>(m.positions[i].y));
}

} // namespace

MinimalGraphSolver::MinimalGraphSolver(const Mesh &m,
			       std::function<double(const Vec2d &)> func)
	: m(m) // mesh reference
	, N(m.vertex_count()) // number of vertices
	, N_b(m.boundary_count()) // number of boundary vertices
	, u(N) // solution vector at current iteration
	, uold(N) // solution vector at previous iteration
	, du(N) // update of solution vector
	, q(m.triangle_count()) // denominator of the nonlinearity at each triangle
	, f(N) // right hand side vector
	, b(N, 0.0) // right hand side vector modified to account for boundary conditions
	, r(N)  // residual vector for CG solver
	, p(N) // search direction vector for CG solver
	, Ap(N) // temporary vector for CG solver
	, inited(false)
	, iterate_N(0)
	, iterate_P(0)
	, residual_N(iter_max, 0.0)
	, residual_P(iter_max, 0.0)
	, converged(false)
{
	// pattern for stiffness matrix of P1 elements
	// only upper triangular part including diagonal since the matrix is symmetric
	// size(P.row_start) = number of vertices + 1
	// size(P.col) = number of non-zero entries in the upper triangular part of the stiffness matrix
	build_P1_CSRPattern(m, P); 

	for (size_t i = 0; i < N; ++i) {
		Vec2d pos(static_cast<double>(m.positions[i].x),
			  static_cast<double>(m.positions[i].y));
		// right hand side at every vertex
		f[i] = func(pos); 
	}
}

void MinimalGraphSolver::clear_solution(bool Newton)
{
	memset(b.data, 0, N * sizeof(double)); // set b to zero
	memset(du.data, 0, N * sizeof(double)); // set du to zero

	// modify b to account for boundary conditions: for each boundary vertex bi, set b[bi] = f[bi]
	for (size_t i = 0; i < N_b; ++i) {
		uint32_t bi = m.boundary[i]; 
		b[bi] = f[bi];
	}

	memcpy(u.data, b.data, N * sizeof(double));  // set u to be with the boundary conditions
	memcpy(uold.data, b.data, N * sizeof(double)); // set uold to b as well
	converged = false;
	inited = true;
	if (Newton) {
		iterate_N = 0;
	} else {
		iterate_P = 0;
	}
}

double MinimalGraphSolver::compute_denominator(TArray<double> &den,
				      const TArray<double> &u)
{
	double area = 0.0;
	// compute denomintor of the non linearity at each trinagle and store in den
	// at the same time compute the area of the surface corresponding to the current solution u, which is the sum of the area of each triangle divided by the corresponding denominator
	for (size_t t = 0; t < m.triangle_count(); ++t) {

		// get vertices of the triangle
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];

		Vec2d A = xy(m, a); // position of vertex a
		Vec2d B = xy(m, b); // position of vertex b
		Vec2d C = xy(m, c); // position of vertex c
		Vec2d AB = {B[0] - A[0], B[1] - A[1]}; // vector from A to B
		Vec2d AC = {C[0] - A[0], C[1] - A[1]}; // vector from A to C

		// compute the area of the triangle
		double ABAB = dot(AB, AB);
		double ACAC = dot(AC, AC);
		double ABAC = dot(AB, AC);
		double tri_area = 0.5 * sqrt(ABAB * ACAC - ABAC * ABAC);
		
		// mult = 1 / (4 * area of the triangle)
		double mult = 0.25 / tri_area; 
		ABAB *= mult; 
		ACAC *= mult;
		ABAC *= mult;

		// compute local stiffness matrix for the triangle defined by AB and AC
		double S_loc[6];
		S_loc[0] = ACAC - 2 * ABAC + ABAB;
		S_loc[1] = ACAC;
		S_loc[2] = ABAB;
		S_loc[3] = ABAC - ACAC;
		S_loc[4] = -ABAC;
		S_loc[5] = ABAC - ABAB;

		den[t] = 1.0 / sqrt(1 
			+ u[a] * u[a] * S_loc[0] 
			+ u[b] * u[b] * S_loc[1] 
			+ u[c] * u[c] * S_loc[2] 
			+ 2 * (u[a] * u[b] * S_loc[3] 
			+ u[b] * u[c] * S_loc[4] 
			+ u[c] * u[a] * S_loc[5]));
		area += tri_area / den[t];
	}
	return area;
}

void MinimalGraphSolver::do_iterate_Newton(size_t max_iter, double tol,
				  const double min_alpha, const double c,
				  const double rho)
{
	if (!inited) {
		clear_solution(true);
	}

	TArray<double> u_tmp(N, 0.0);
	int iterCG;
	double error2 = 0.0; // squared norm of the update du
	double errorCG = 0.0; // residual error of the CG solver
	double area; // area of the surface corresponding to the current solution u
	double alpha = 1.0; // step size for the line search, initialized to 1.0
	double energy_tmp = 0.0; 
	bool flag = true;
	double prev_error = 1e30;
	int stagnant_count = 0;

	// q is denominator of the nonlinearity at each triangle, computed based on the current solution u
	area = compute_denominator(q, u);
	if (iterate_N == 0) {
		printf("Starting Newton solver.... \n");
		printf("%-10s %-15s %-15s %-15s %-15s\n", "Iter", "ErrorNewton",
		       "IterCG", "ErrorCG", "Area");
		printf("%-10s %-15s %-15s %-15s %-15g \n", "-", "-", "-", "-",
		       area);
	}

	size_t target_iter = iterate_N + max_iter;
	while (iterate_N < target_iter) {
		build_P1_stiffness_matrix_NS(m, P, S_modified, q.data, u.data, area);
		build_P1_rhs_NS(m, q.data, u.data, b);

		for (size_t i = 0; i < N_b; ++i) {
			uint32_t bi = m.boundary[i];
			S_modified(bi, bi) = HUGE_DIAG;
			b[bi] = 0;
		}

		iterCG = conjugate_gradient_solve(S_modified, b.data, du.data, r.data,
						 p.data, Ap.data, &errorCG,
						 tolCG, 10000, false);

		area = compute_denominator(q, u);
		error2 = blas_dot(du.data, du.data, N);
		flag = true;
		alpha = 1.0;

		while (flag) {
			for (size_t i = 0; i < N; ++i) {
				u_tmp[i] = u[i] + du[i];
			}
			energy_tmp = compute_denominator(q, u_tmp);
			if (energy_tmp < area - c * alpha * error2 || alpha <= min_alpha) {
				flag = false;
			} else {
				alpha *= rho;
			}
		}

		for (size_t i = 0; i < N; ++i) {
			u[i] += alpha * du[i];
		}

		error2 = sqrt(error2);
		printf("%-10ld %-15g %-15d %-15g %-15g\n", iterate_N, error2,
		       iterCG, errorCG, area);

		if (iterate_N < residual_N.size) {
			residual_N[iterate_N] = error2;
		}
		++iterate_N;

		if (error2 < tol) {
			converged = true;
			printf("Converged after %ld iterations.\n", iterate_N);
			break;
		}

		if (error2 > prev_error * 0.9999) {
			++stagnant_count;
			if (stagnant_count > 10) {
				printf("Stagnation detected after %ld iterations (error not decreasing).\n",
				       iterate_N);
				converged = false;
				break;
			}
		} else {
			stagnant_count = 0;
		}
		prev_error = error2;

		if (iterate_N >= iter_max) {
			printf("Reached maximum iteration limit (%zu).\n", iter_max);
			break;
		}
	}

	if (!converged && error2 <= tol) {
		converged = true;
		printf("Converged after %ld iterations.\n", iterate_N);
	} else if (!converged) {
		printf("Did not converge after %ld iterations.\n", iterate_N);
	}
}

void MinimalGraphSolver::do_iterate_Picardi(size_t max_iter, double tol)
{
	(void)max_iter;
	(void)tol;
	/* Initial skeleton: Picard iteration is intentionally left as a follow-up. */
}
