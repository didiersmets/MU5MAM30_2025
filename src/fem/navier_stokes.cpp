#include "fem/navier_stokes.h"

#include "common/vec3.h"
#include "fem/P1.h"
#include "linalg/conjugate_gradient.h"
#include "linalg/tiny_blas.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <vector>

NavierStokesSolver::NavierStokesSolver(const Mesh& m)
    : m(m), N(m.vertex_count()), omega(N), Momega(N), psi(N), r(N), p(N), Ap(N)
{
  build_P1_CSRPattern(m, P);
  build_P1_mass_matrix(m, P, M);
  build_P1_stiffness_matrix(m, P, S);

  vol    = M.sum();  // volume
  inited = false;
  t      = 0;

  /// current residue r = Mf - Su
}

// we need to be sure that the solution PSI coming from The poisson solver S_PSI = - M_OMEGA
// has zero mean value over the domain
// that means that the integrale over the domain of PSI is zero

void NavierStokesSolver::set_zero_mean(double* V)
{
  /* Your implementation goes here */

  // In order to compute the integral over the domain of V, we use the mass matrix M
  // as int(V) = sum i sum j M_ij * V_j
  // recall that Sum i (phi i) = 1
  // integral = sum i (Vi * int(phi_i)) = sum i (Vi * int( 1 * phi_i))=sum i (Vi * sum j M_ij) = sum
  // i (M * V)_i

  double integral = 0.0;
  M.mvp(V, Ap.data);  // Momega = M * V

  for (size_t i = 0; i < N; i++)
  {
    integral += Ap.data[i];
  }

  double mean_value = integral / vol;

  // we subtract the mean value to each entry of V to ensure zero mean
  for (size_t i = 0; i < N; i++)
  {
    V[i] -= mean_value;
  }
}

void compute_coriolis(const Mesh& m, double* coriolis, double omega_earth)
{
  // const double omega_earth = 2.0 / M_PI;  // angular velocity of the earth
  size_t N = m.vertex_count();
  for (size_t i = 0; i < N; i++)
  {
    // On a unit sphere, sin(latitude) corresponds to the z-coordinate
    double z    = m.positions[i].z;
    coriolis[i] = 2.0 * omega_earth * z;
  }
}

/* To compute the TRANSPORT TERM We use the formula :
 \forall j \in I, T[j] = \sum_{i, k} \Omega_i * \Psi_k \int_{\Omega} \phi_i * (\nabla^T \phi_k .
 \nabla \phi_j)

I used the change of variable to the reference triangle (0,0) (1,0) (0,1)
the integral on the reference triangle of phi_i = 1/6 (notice that gradient terms are constant
vector for P1 elements) det(J) = 2 * Area(T_ABC)
*/
void NavierStokesSolver::compute_transport(double* T)
{
  memset(T, 0, N * sizeof(double));  // transport_term

  double* OMEGA = omega.data;
  double* PSI   = psi.data;

  size_t nt = m.triangle_count();

  /* On each triangle we have the contribution to T[i] for the three vertices of the triangle
    1) The integral over the triangle of phi_i = Area / 6
    2) The term (grad_purp \phi_k . \grad \phi_j) is constant over the triangle since P1 elements
   are
   * The gradient of the basis functions \nabla \phi_i are of order 1/L (L = length scale of the
   * triangle) term of order 1/Area (since 1/L * 1/L = 1/L^2 ~ 2/Area). The Area terms cancel out
   * perfectly, leaving a purely topological constant (1/6). Differently from what I had previously
   computed*/

  for (size_t tri = 0; tri < nt; tri++)
  {
    uint32_t a = m.indices[3 * tri];
    uint32_t b = m.indices[3 * tri + 1];
    uint32_t c = m.indices[3 * tri + 2];

    /*
=========================================================================
       PROJECT UPDATE: CORIOLIS FORCE IMPLEMENTATION
       =========================================================================
       We modify the transport term to conserve Absolute Vorticity (eta) instead
       of just Relative Vorticity (omega).

       1. Transport Equation:
          The term J(psi, omega) becomes J(psi, omega + f).

       2. Coriolis Parameter 'f':
          Formula: f = 2 * Omega_earth * sin(latitude)
          Geometry: On a Unit Sphere (Radius = 1), sin(latitude) corresponds exactly
                    to the z-coordinate.
          Therefore: f = 2 * Omega_earth * z

          Instead of computing two separate integrals (one for omega, one for f),
          we can simply sum the scalar values at the nodes first:
          eta_node = omega_node + f_node

        Because the FEM integration logic T(.) is linear with respect to the
        coefficients, passing this sum is mathematically equivalent:
          T(omega + f) == T(omega) + T(f)

       3. Scaling (Omega_earth):
          The professor confirmed that for a unit sphere, the rotation rate is a
          tunable scaling parameter.
          - Test with Omega_earth = 0.0 for standard isotropic turbulence.
          - Test with Omega_earth ~ 10.0 - 50.0 to observe Rossby waves and zonal jets.
      */
    assert(a < N && b < N && c < N);
    double omega_sum = OMEGA[a] + OMEGA[b] + OMEGA[c];
    T[a] += (omega_sum * (PSI[b] - PSI[c])) / 6.0;
    T[b] += (omega_sum * (PSI[c] - PSI[a])) / 6.0;
    T[c] += (omega_sum * (PSI[a] - PSI[b])) / 6.0;
  }
}

// To compute the stream function PSI from the vorticity OMEGA we need to solve the linear system
// associated to the poisson problem S * PSI = - M * OMEGA
size_t NavierStokesSolver::compute_stream_function()
{
  size_t iter = 0;

  // we first compute M * OMEGA
  M.mvp(omega.data, Momega.data);  // Momega = M * omega

  // we set the right hand side b = - M * OMEGA
  for (size_t i = 0; i < N; i++)
  {
    Momega.data[i] = -Momega.data[i];
  }

  double rel_error      = 0.0;
  double tol            = 1e-6;
  int    max_iterations = 1000;
  // solve with conjugate gradient S * PSI = b
  iter = conjugate_gradient_solve(S,
                                  Momega.data,
                                  psi.data,
                                  r.data,
                                  p.data,
                                  Ap.data,
                                  &rel_error,
                                  tol,
                                  max_iterations);

  return iter;
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
  // compute_transport_coriolis(transport.data, double omega_earth);

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
  A.rows      = M.rows;
  A.cols      = M.cols;
  A.nnz       = M.nnz;
  A.symmetric = M.symmetric;
  A.row_start = M.row_start;
  A.col       = M.col;
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

  double rel_error      = 0.0;
  double tolerance      = 1e-6;
  int    max_iterations = 1000;

  size_t iterations = conjugate_gradient_solve(A,
                                               Momega.data,
                                               omega.data,
                                               r.data,
                                               p.data,
                                               Ap.data,
                                               &rel_error,
                                               tolerance,
                                               max_iterations,
                                               false);

  // Necessary to ensure that omega has zero mean at each time step
  set_zero_mean(omega.data);

  t += dt;
}

void NavierStokesSolver::compute_transport_coriolis(double* T, double omega_earth)
{
  memset(T, 0, N * sizeof(double));  // transport_term
  std::vector<double> coriolis(N);
  compute_coriolis(m, coriolis.data(), omega_earth);

  double* OMEGA = omega.data;
  double* PSI   = psi.data;

  size_t nt = m.triangle_count();

  /* On each triangle we have the contribution to T[i] for the three vertices of the triangle
    1) The integral over the triangle of phi_i = Area / 6
    2) The term (grad_purp \phi_k . \grad \phi_j) is constant over the triangle since P1 elements
   are
   * The gradient of the basis functions \nabla \phi_i are of order 1/L (L = length scale of the
   * triangle) term of order 1/Area (since 1/L * 1/L = 1/L^2 ~ 2/Area). The Area terms cancel out
   * perfectly, leaving a purely topological constant (1/6). Differently from what I had previously
   computed*/

  for (size_t tri = 0; tri < nt; tri++)
  {
    uint32_t a = m.indices[3 * tri];
    uint32_t b = m.indices[3 * tri + 1];
    uint32_t c = m.indices[3 * tri + 2];

    /*
=========================================================================
       PROJECT UPDATE: CORIOLIS FORCE IMPLEMENTATION
       =========================================================================
       We modify the transport term to conserve Absolute Vorticity (eta) instead
       of just Relative Vorticity (omega).

       1. Transport Equation:
          The term J(psi, omega) becomes J(psi, omega + f).

       2. Coriolis Parameter 'f':
          Formula: f = 2 * Omega_earth * sin(latitude)
          Geometry: On a Unit Sphere (Radius = 1), sin(latitude) corresponds exactly
                    to the z-coordinate.
          Therefore: f = 2 * Omega_earth * z

          Instead of computing two separate integrals (one for omega, one for f),
          we can simply sum the scalar values at the nodes first:
          eta_node = omega_node + f_node

        Because the FEM integration logic T(.) is linear with respect to the
        coefficients, passing this sum is mathematically equivalent:
          T(omega + f) == T(omega) + T(f)

       3. Scaling (Omega_earth):
          The professor confirmed that for a unit sphere, the rotation rate is a
          tunable scaling parameter.
          - Test with Omega_earth = 0.0 for standard isotropic turbulence.
          - Test with Omega_earth ~ 10.0 - 50.0 to observe Rossby waves and zonal jets.
      */
    assert(a < N && b < N && c < N);
    double omega_sum = OMEGA[a] + OMEGA[b] + OMEGA[c];
    omega_sum += coriolis[a] + coriolis[b] + coriolis[c];  // linear addition of coriolis term
    T[a] += (omega_sum * (PSI[b] - PSI[c])) / 6.0;
    T[b] += (omega_sum * (PSI[c] - PSI[a])) / 6.0;
    T[c] += (omega_sum * (PSI[a] - PSI[b])) / 6.0;
  }
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
  compute_transport_coriolis(transport.data, omega_earth);  // transport = T(OMEGA + f, PSI)

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
  A.rows      = M.rows;
  A.cols      = M.cols;
  A.nnz       = M.nnz;
  A.symmetric = M.symmetric;
  A.row_start = M.row_start;
  A.col       = M.col;
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

  /* Your implementation goes here */

  double rel_error      = 0.0;
  double tolerance      = 1e-6;
  int    max_iterations = 1000;

  size_t iterations = conjugate_gradient_solve(A,
                                               Momega.data,  // RHS
                                               omega.data,   // Initial guess AND result
                                               r.data,
                                               p.data,
                                               Ap.data,
                                               &rel_error,
                                               tolerance,
                                               max_iterations,
                                               false);

  // Necessary to ensure that omega has zero mean at each time step
  set_zero_mean(omega.data);

  t += dt;
}