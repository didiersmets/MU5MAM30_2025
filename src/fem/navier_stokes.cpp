#include "navier_stokes.h"

#include "vec3.h"
#include "P1.h"
#include "conjugate_gradient.h"
#include "tiny_blas.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <vector>

NavierStokesSolver::NavierStokesSolver(const Mesh& m)
    : m(m), N(m.vertex_count()), omega(N), Momega(N), psi(N), r(N), p(N), Ap(N), velocity(0)
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

// For the project I introducee the coriolis term 
// so create a function to compute the coriolis term

void compute_coriolis(const Mesh& m, double* coriolis, double omega_earth)
{
  // double omega_earth = 2.0 / M_PI;  // angular velocity of the earth
  size_t N = m.vertex_count();
  for (size_t i = 0; i < N; i++)
  {
    // On a unit sphere, sin(latitude) corresponds to the z-coordinate
    double z    = m.positions[i].z;
    coriolis[i] = 2.0 * omega_earth * z;
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


	// compute transport term contribution for each vertex of the triangle
	// I think that the area term is not necessary since it is included in the mass matrix
    double omega_sum = OMEGA[a] + OMEGA[b] + OMEGA[c];
    T[a] += (omega_sum * (PSI[c] - PSI[b])) / 6.0 ;
    T[b] += (omega_sum * (PSI[a] - PSI[c])) / 6.0 ;
    T[c] += (omega_sum * (PSI[b] - PSI[a])) / 6.0 ;
  }
}

// function to use to change the variables of the initial condition on the fluid and see better the coriolis effect
// 100 * x * exp(-50*x^2) * (1 + 0.5 * cos(0.05 * atan2(z, y)))

// Now we compute the transport term with the coriolis effect included
void NavierStokesSolver::compute_transport_coriolis(double* T , double omega_earth)
{
  memset(T, 0, N * sizeof(double));              // transport_term
   // angular velocity of the earth
  std::vector<double> coriolis(N);
  compute_coriolis(m, coriolis.data(), omega_earth);

  double* OMEGA = omega.data;
  double* PSI   = psi.data;

  size_t nt = m.triangle_count();


  for (size_t tri = 0; tri < nt; tri++)
  {
    uint32_t a = m.indices[3 * tri];
    uint32_t b = m.indices[3 * tri + 1];
    uint32_t c = m.indices[3 * tri + 2];


    assert(a < N && b < N && c < N);
    double omega_sum = (OMEGA[a] + OMEGA[b] + OMEGA[c]) / 6.0;      // here we have the integral over the reference triangle of phi_i which is 1/6
    omega_sum += (coriolis[a] + coriolis[b] + coriolis[c]) / 6.0;   // here we have the integral over the reference triangle of 1 which is 1/2
    T[a] += (omega_sum * (PSI[b] - PSI[c]));
    T[b] += (omega_sum * (PSI[c] - PSI[a]));
    T[c] += (omega_sum * (PSI[a] - PSI[b]));
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
  //compute_transport_coriolis(transport.data, double omega_earth);

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

  [[maybe_unused]] size_t iterations = conjugate_gradient_solve(A, Momega.data, omega.data, r.data, p.data, Ap.data, &rel_error, tolerance, max_iterations, false);

  // Necessary to ensure that omega has zero mean at each time step
  set_zero_mean(omega.data);


  t += dt;
}

void NavierStokesSolver::time_step_coriolis(double dt, double nu, double omega_earth)
{

  // first computation of PSI comes from an omega with zero mean: see src/test_navier_stokes.cpp
  // compute PSI from OMEGA
  compute_stream_function();

  //  we ensure that PSI over the domain is zero too.
  set_zero_mean(psi.data);

  // compute transport term T(OMEGA, PSI)
  TArray<double> transport(N);
  compute_transport_coriolis(transport.data, omega_earth);

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

  [[maybe_unused]] size_t iterations = conjugate_gradient_solve(A, Momega.data, omega.data, r.data, p.data, Ap.data, &rel_error, tolerance, max_iterations, false);

  // Necessary to ensure that omega has zero mean at each time step
  set_zero_mean(omega.data);


  t += dt;
}