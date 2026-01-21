#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"
#include "conjugate_gradient.h"

#include "P1.h"
#include "tiny_blas.h"

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
	, r(N)
	, p(N)
	, Ap(N)
  , T(N)
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);

  A.rows      = M.rows;
  A.cols      = M.cols;
  A.nnz       = M.nnz;
  A.symmetric = M.symmetric;
  A.row_start = M.row_start;
  A.col       = M.col;
  A.data.resize(A.nnz);
#endif
	vol = M.sum();
	inited = false;
	t = 0;
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	/* Your implementation goes here */
  double integral = 0.0;
  TArray<double> MV(N);
  M.mvp(V, MV.data);
  for (size_t i = 0; i < N; i++) {
      integral += MV[i];
  }
  const double mean = integral / vol;
  for (size_t i = 0; i < N; i++) {
      V[i] -= mean;
  }

}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */
	double* v_omega = omega.data;
	double* v_psi = psi.data;

	/*

	T(Omega(t), Psi(t))_k = Sum_{i,j \in [1,N]} Omega_i(t) * Psi_j(t) * \integral psi_i * grad_perp(psi_j) \dot grad(psi_k) dA

	The two gradients are constants and the integral of psi is 1/6.
	T_A = 2 * A(T) / 6 * (\sum_{i in triangle} Omega_i) * (\sum_{j in triangle} Psi_j * grad_perp(psi_j) \dot grad(psi_A))
						 ^-- is computed once [sum_omega]						^						   ^-- is Zero for j = A, +1 or -1 in the other cases
																				^-- is stored in the v_psi array
	
																				The integral over the whole domain can be splitted as a sum of integrals over each triangle.
	We conpute this on each triangle and accumulate the results.

	The computation is done on a standard triangle, the transofrmation is up to a constant factor.

	T_i =
	i \in {A,B,C}
	
	
	*/

	// iterate over all triangles
	size_t tri_count = m.triangle_count();
	for (size_t t = 0; t < tri_count; t++) {
		// get vertex indices
		uint32_t ia = m.indices[3 * t + 0];
		uint32_t ib = m.indices[3 * t + 1];
		uint32_t ic = m.indices[3 * t + 2];

		// get vertex positions
		Vec3f A = m.positions[ia];
		Vec3f B = m.positions[ib];
		Vec3f C = m.positions[ic];

		// compute edges
		Vec3f AB = B - A;
		Vec3f AC = C - A;

		double sum_omega = v_omega[ia] + v_omega[ib] + v_omega[ic];
		// compute local contributions to transport term
		T[ia] += sum_omega * (v_psi[ic] - v_psi[ib]) / 6.0f;
		T[ib] += sum_omega * (v_psi[ia] - v_psi[ic]) / 6.0f;
		T[ic] += sum_omega * (v_psi[ib] - v_psi[ia]) / 6.0f;
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/* Your implementation goes here */

  // compute Momega
  M.mvp(omega.data, Momega.data);

  // rhs is -Momega
  for (size_t i = 0; i < N; i++)
  {
    Momega.data[i] = -Momega.data[i];
  }

  double rel_error = 0.0;
  iter = conjugate_gradient_solve(S,
                                  Momega.data,
                                  psi.data,
                                  r.data,
                                  p.data,
                                  Ap.data,
                                  &rel_error,
                                  tol,
                                  iter_max);

  // restore value of Momega
  for (size_t i = 0; i < N; i++)
  {
    Momega.data[i] = -Momega.data[i];
  } 

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function();
  set_zero_mean(psi.data);

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Your implementation goes here */

  // compute transport term
  compute_transport(T.data);

  // compute rhs = M * omega + dt * T
  TArray<double> rhs(N);
  for(size_t i = 0; i < N; i++) {
      rhs[i] = Momega[i] + dt * T[i];
  }

  // compute matrix A = M + nu * dt * S (only once)
  if(!inited) {
    for(size_t i = 0; i < A.nnz; i++) {
        A.data[i] = M.data[i] + nu * dt * S.data[i];
    }
    inited = true;
  }

  // solve system
  double rel_error = 0.0;
  size_t iterations = conjugate_gradient_solve(A,
                                               rhs.data,
                                               omega.data,
                                               r.data,
                                               p.data,
                                               Ap.data,
                                               &rel_error,
                                               tol,
                                               iter_max,
                                               false);



	set_zero_mean(omega.data);

	t += dt;
}
