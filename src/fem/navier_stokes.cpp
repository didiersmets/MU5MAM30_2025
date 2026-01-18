#include "navier_stokes.h"

#include "vec3.h"
#include "P1.h"
#include "conjugate_gradient.h"
#include "tiny_blas.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>

NavierStokesSolver::NavierStokesSolver(const Mesh& m)
    : m(m), N(m.vertex_count()), omega(N), Momega(N), psi(N), r(N), p(N), Ap(N)
{
  build_P1_CSRPattern(m, P);
  build_P1_mass_matrix(m, P, M);
  build_P1_stiffness_matrix(m, P, S);

  vol = M.sum();  // volume
  inited = false;
  t = 0;

}


void NavierStokesSolver::set_zero_mean(double* V)
{
  /* Your implementation goes here */

  double integral = 0.0;
  M.mvp(V, Momega.data);  // Momega = M * V

  for (size_t i = 0; i < N; i++)
  {
    integral += Momega.data[i];
  }

  double mean_value = integral / vol;

  // we subtract the mean value to each entry of V to ensure zero mean
  for (size_t i = 0; i < N; i++)
  {
    V[i] -= mean_value;
  }
}


void NavierStokesSolver::compute_transport(double* T)
{
  memset(T, 0, N * sizeof(double));  // transport_term

  double* OMEGA = omega.data;
  double* PSI   = psi.data;

  size_t nt = m.triangle_count();

  for (size_t tri = 0; tri < nt; tri++)
  {
    uint32_t a = m.indices[3 * tri];
    uint32_t b = m.indices[3 * tri + 1];
    uint32_t c = m.indices[3 * tri + 2];

    double omega_sum = OMEGA[a] + OMEGA[b] + OMEGA[c];
    T[a] += (omega_sum * (PSI[c] - PSI[b])) / 6.0;
    T[b] += (omega_sum * (PSI[a] - PSI[c])) / 6.0;
    T[c] += (omega_sum * (PSI[b] - PSI[a])) / 6.0;
  }
}

size_t NavierStokesSolver::compute_stream_function()
{
  size_t iteration = 0;

  /* Your implementation goes here */
  // we first compute M * OMEGA
  M.mvp(omega.data, Momega.data);  // Momega = M * omega

  // we set the right hand side b = - M * OMEGA
  for (size_t i = 0; i < N; i++)
  {
    Momega.data[i] = -Momega.data[i];
  }

  double relative_error = 0.0;
  double tolerance = 1e-6;
  int    max_iterations = 1000;
  // solve with conjugate gradient S * PSI = b

  iteration = conjugate_gradient_solve(S,Momega.data, psi.data, r.data, p.data, Ap.data, &relative_error, tolerance, max_iterations);

  return iteration;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
  // first computation of PSI comes from an omega with zero mean: see src/test_navier_stokes.cpp
  // compute PSI from OMEGA
  compute_stream_function();
  //  we ensure that PSI over the domain is zero too.
  set_zero_mean(psi.data);

  // compute transport term T(OMEGA, PSI)
  TArray<double> transport(N);
  compute_transport(transport.data);

  // compute the right hand side
  // rhs = M * OMEGA + dt * T(OMEGA, PSI)
  M.mvp(omega.data, Momega.data);  // Momega = M * omega

  // Momega = dt * transport + Momega
  // tiny_blas.axpy(N, dt, transport.data, Momega.data);
  for (size_t i = 0; i < N; i++)
  {
    Momega.data[i] += dt * transport.data[i];
  }

  // build the left hand side matrix A = M + nu * dt * S
  CSRMatrix A;
  A.rows = M.rows;
  A.cols = M.cols;
  A.nnz = M.nnz;
  A.symmetric = M.symmetric;
  A.row_start = M.row_start;
  A.col = M.col;
  A.data.resize(A.nnz);

  // compute A = M + S
  double factor = nu * dt;
  for (size_t k = 0; k < A.nnz; k++)
  {
    A.data[k] = M.data[k] + ((factor) *S.data[k]);
  }

  /**********************************************************************
   * Solve the system :
   *
   *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
   *
   *********************************************************************/


  double rel_error = 0.0;
  double tolerance = 1e-6;
  int    max_iterations = 1000;

  size_t iterations = conjugate_gradient_solve(A, Momega.data, omega.data, r.data, p.data, Ap.data, &rel_error, tolerance, max_iterations, false);

  // Necessary to ensure that omega has zero mean at each time step
  set_zero_mean(omega.data);


  t += dt;
}