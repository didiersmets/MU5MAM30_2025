#include "fem/navier_stokes.h"

#include "fem/P1.h"
#include "linalg/conjugate_gradient.h"
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

  /* Your implementation goes here */

  double* OMEGA = omega.data;
  double* PSI   = psi.data;

  size_t nt = m.triangle_count();

  /*
  I used the change of variable to the reference triangle (0,0) (1,0) (0,1)
  the integral on the triangle of phi_i = 1/6
  det(J) = 2* Area(T_ABC)
  */

  /* We use the formula :
   \forall j \in I, T[j] = \sum_{i, k} \Omega_i * \Psi_k \int_{\Omega} \phi_i * (\nabla^T \phi_k .
   \nabla \phi_j) */

  for (size_t tri = 0; tri < nt; tri++)
  {
    uint32_t a = m.indices[3 * tri];
    uint32_t b = m.indices[3 * tri + 1];
    uint32_t c = m.indices[3 * tri + 2];

    Vec3f A_pos = m.positions[a];
    Vec3f B_pos = m.positions[b];
    Vec3f C_pos = m.positions[c];
    // Edge vectors
    Vec3d AB = {(double) B_pos[0] - A_pos[0],
                (double) B_pos[1] - A_pos[1],
                (double) B_pos[2] - A_pos[2]};
    Vec3d AC = {(double) C_pos[0] - A_pos[0],
                (double) C_pos[1] - A_pos[1],
                (double) C_pos[2] - A_pos[2]};

    // Calculate area of the triangle
    double area = 0.5 * norm(cross(AB, AC));

    double omega_sum = OMEGA[a] + OMEGA[b] + OMEGA[c];
    T[a] += (area / 3.0) * omega_sum * (PSI[c] - PSI[b]);
    T[b] += (area / 3.0) * omega_sum * (PSI[a] - PSI[c]);
    T[c] += (area / 3.0) * omega_sum * (PSI[b] - PSI[a]);
  }
}

// To compute the stream function PSI from the vorticity OMEGA we need to solve the linear system
// associated to the poisson problem S * PSI = - M * OMEGA

size_t NavierStokesSolver::compute_stream_function()
{
  size_t iter = 0;

  /* Your implementation goes here */
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
  compute_stream_function();

  /**********************************************************************
   * Solve the system :
   *
   *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
   *
   *********************************************************************/

  /* Your implementation goes here */

  // we ensure that the vorticity omega has zero mean at each time step
  // hence, we are sure that integral of
  //  PSI over the domain is zero too.
  set_zero_mean(omega.data);

  t += dt;
}