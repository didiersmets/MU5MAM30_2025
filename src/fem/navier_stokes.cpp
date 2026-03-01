#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "conjugate_gradient.h"
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
  , RHS(N)
  , other(N)
{
#if USE_FEM_MATRIX
  build_P1_mass_matrix(m, M);
  build_P1_stiffness_matrix(m, S);
#else
  build_P1_CSRPattern(m, P);
  build_P1_mass_matrix(m, P, M);
  build_P1_stiffness_matrix(m, P, S);
#endif
  vol = M.sum();
  inited = false;
  t = 0;
}

void NavierStokesSolver::set_zero_mean(double *V)
{
  /* (done) Your implementation goes here */
  M.mvp(V, other.data);
  double c = blas_sum_in_place(other.data, N) / vol;
  for (size_t i=0; i<N; ++i)
    V[i] -= c;
}

void NavierStokesSolver::compute_transport(double *T)
{
  memset(T, 0, N * sizeof(double));

  /* (done) Your implementation goes here */
  size_t n_tri = m.triangle_count();

  for (size_t tri=0; tri<n_tri; ++tri) {
    uint32_t va = m.indices[3*tri];
    uint32_t vb = m.indices[3*tri+1];
    uint32_t vc = m.indices[3*tri+2];

    double omega_sum = omega[va] + omega[vb] + omega[vc];

    T[va] += (psi[vc] - psi[vb]) * omega_sum / 6.0;
    T[vb] += (psi[va] - psi[vc]) * omega_sum / 6.0;
    T[vc] += (psi[vb] - psi[va]) * omega_sum / 6.0;
  }
}

size_t NavierStokesSolver::compute_stream_function()
{
  size_t iter = 0;

  /* (done) Your implementation goes here */
  M.mvp(omega.data, Momega.data);
  double rel_error;

  iter = conjugate_gradient_solve(S,
				  Momega.data,
				  psi.data,
				  r.data,
				  p.data,
				  Ap.data,
				  &rel_error,
				  tol,
				  iter_max,
				  inited);
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

  /* (done) Your implementation goes here */
  double *omegad = omega.data;
  double *Momegad = Momega.data;
  double *rd = r.data;
  double *pd = p.data;
  double *Apd = Ap.data;
  double *bd = RHS.data;
  double *otherd = other.data;

  double rel_error, r2, b2;

  size_t iter_cg = 0;

  /* We need to rewrite the conjugate gradient algorithm here in order to not have to
     create the matrix M + \nu * dt * S each time */

  /* Initialize CG variables */
  /* b = M * omega(t) + dt * T(Omega,Psi)(t) */
  compute_transport(bd);
  blas_axpby(1, Momegad, dt, bd, N);
  b2 = blas_dot(bd, bd, N);
  /* r_0 = b - Ax_0 */
  S.mvp(omegad, rd);
  blas_axpby(1, Momegad, nu*dt, rd, N);
  blas_axpby(1, bd, -1, rd, N);
  r2 = blas_dot(rd, rd, N);
  /* p_0 = r_0 */
  blas_copy(rd, pd, N);

  rel_error = sqrt(r2 / b2);


  while ((iter_cg < iter_max) && (rel_error > tol)) {
    /* alpha_k = r_k^Tr_k / (p_k^T A p_k) */
    M.mvp(pd, Apd);
    S.mvp(pd, otherd);
    blas_axpby(nu*dt, otherd, 1, Apd, N);
    double alpha = r2 / blas_dot(pd, Apd, N);

    /* x_{k+1} = x_k + \alpha_k p_k */
    blas_axpy(alpha, pd, omegad, N);
    /* r_{k+1} = r_k - \alpha_k Ap_k*/
    blas_axpy(-alpha, Apd, rd, N);

    /* r2_new = r_{k+1}^T r_{k+1} */
    double r2_new = blas_dot(rd, rd, N);

    /* beta_k = r_{k+1}^T r_{k+1} / (r_k^T r_k) */
    /* p_{k+1} = r_{k+1} + beta_{k+1} p_k */
    double beta = r2_new / r2;
    blas_axpby(1, rd, beta, pd, N);

    r2 = r2_new;
    rel_error = sqrt(r2 / b2);

    M.mvp(omegad, Momegad);

    iter_cg++;
  }

  set_zero_mean(omega.data);

  t += dt;
}
