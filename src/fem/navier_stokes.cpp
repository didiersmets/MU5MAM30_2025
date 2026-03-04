#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"
#include "conjugate_gradient.h"
#include "tiny_blas.h"

#if USE_P2
  #include "triangle_quadrature.h"
  #include "P2.h"
  #include <algorithm>
#else
  #include "P1.h"
#endif

static size_t inline dof_count(const Mesh &m) {
#if USE_P2
  return  m.vertex_count() + m.index_count() / 2;
#else
  return  m.vertex_count();
 #endif
}

#if USE_P2
static auto get_transport_integrand_P2(const NavierStokesSolver *NS,
				     const uint32_t dof[6],
				     const Vec3 grad_lambda[3],
				     const Vec3 *p_n,
				     float phi[6],
				     Vec3 grad_phi[6],
				     Vec3 *p_grad_psi_perp,
				     uint32_t k)
{
  /* Return the quadrature integrand function at barycentric coordinates (la,lb,lc) */
  const double *psi = NS->psi.data;
  const double *omega = NS->omega.data;

  auto integrand = [=](float la, float lb, float lc) {
    /* Get the P2 shape functions */
    phi[0] = la*(2*la-1);
    phi[1] = lb*(2*lb-1);
    phi[2] = lc*(2*lc-1);
    phi[3] = 4*la*lb;
    phi[4] = 4*lb*lc;
    phi[5] = 4*lc*la;

    /* Get the P2 shape function gradients */
    grad_phi[0].dcopy((4*la-1)*grad_lambda[0]);
    grad_phi[1].dcopy((4*lb-1)*grad_lambda[1]);
    grad_phi[2].dcopy((4*lc-1)*grad_lambda[2]);
    grad_phi[3].dcopy(float(4)*(la*grad_lambda[1] + lb*grad_lambda[0]));
    grad_phi[4].dcopy(float(4)*(lb*grad_lambda[2] + lc*grad_lambda[1]));
    grad_phi[5].dcopy(float(4)*(lc*grad_lambda[0] + la*grad_lambda[2]));

    /* Get grad psi perp */
    p_grad_psi_perp->dcopy(Vec3::Zero);

    for (uint32_t i=0; i<6; ++i)
      *p_grad_psi_perp += float(psi[dof[i]]) * grad_phi[i];

    *p_grad_psi_perp = cross(*p_n, *p_grad_psi_perp);
    double dot_p = dot(*p_grad_psi_perp, grad_phi[k]);

    /* Get the value of omega */
    double omega_val = 0;
    for (uint32_t i=0; i<6; ++i)
      omega_val += omega[dof[i]] * phi[i];

    return omega_val * dot_p;
  };
  return integrand;
}


static void compute_transport_P2(const NavierStokesSolver *NS, double *T)
{
  const Mesh &m = NS->m;
  const size_t n_tri = m.triangle_count();

  for (size_t tri = 0; tri < n_tri; ++tri) {
    /* Get the dofs */
    uint32_t v_ids[3] = {m.indices[3*tri],
			 m.indices[3*tri+1],
			 m.indices[3*tri+2]};
    std::sort(std::begin(v_ids), std::end(v_ids));
    uint32_t va = v_ids[0];
    uint32_t vb = v_ids[1];
    uint32_t vc = v_ids[2];
    uint32_t eab = m.edge2dof.find(pack(va, vb))->second;
    uint32_t eac = m.edge2dof.find(pack(va, vc))->second;
    uint32_t ebc = m.edge2dof.find(pack(vb, vc))->second;
    /* Careful: different order in dofs ! */
    uint32_t dof[6] = {va, vb, vc, eab, ebc, eac};

    /* Get the geometry */
    Vec3f A = m.positions[va];
    Vec3f B = m.positions[vb];
    Vec3f C = m.positions[vc];

    Vec3 AB = B - A;
    Vec3 AC = C - A;
    Vec3 Nvec = cross(AB,AC);
    double norm_n = norm(Nvec);
    Vec3 n = Nvec / float(norm_n);

    /*  Compute the gradients of the barycentric coordinates */
    Vec3 grad_lambda[3];
    grad_lambda[0] = cross(n, C - B) / float(norm_n);
    grad_lambda[1] = cross(n, A - C) / float(norm_n);
    grad_lambda[2] = cross(n, B - A) / float(norm_n);

    /* Fill T */
    float phi[6];
    Vec3 grad_phi[6];
    Vec3 grad_psi_perp;

    for (uint32_t k=0; k<6; ++k) {
      auto integrand = get_transport_integrand_P2(NS,
						  dof,
						  grad_lambda,
						  &n,
						  phi,
						  grad_phi,
						  &grad_psi_perp,
						  k);
      T[dof[k]] += integrate_triangle_deg5(A, B, C, integrand);
    }
  }
}
#else
static void compute_transport_P1(const NavierStokesSolver *NS, double *T)
{
  size_t n_tri = NS->m.triangle_count();

  for (size_t tri=0; tri<n_tri; ++tri) {
    uint32_t va = NS->m.indices[3*tri];
    uint32_t vb = NS->m.indices[3*tri+1];
    uint32_t vc = NS->m.indices[3*tri+2];

    double omega_sum = NS->omega[va] + NS->omega[vb] + NS->omega[vc];

    T[va] += (NS->psi[vc] - NS->psi[vb]) * omega_sum / 6.0;
    T[vb] += (NS->psi[va] - NS->psi[vc]) * omega_sum / 6.0;
    T[vc] += (NS->psi[vb] - NS->psi[va]) * omega_sum / 6.0;
  }
}
#endif

NavierStokesSolver::NavierStokesSolver(Mesh &m)
  : m(m)
  , N(dof_count(m))
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
#elif USE_P2
  build_P2_CSRPattern(m, P);
  build_P2_mass_matrix(m, P, M);
  build_P2_stiffness_matrix(m, P, S);

  uint32_t n_vtx = m.vertex_count();
  uint32_t n_edges =  m.index_count() / 2;
  m.e2vtx.resize(2 * n_edges);
  for (const auto &elt : m.edge2dof) {
    uint32_t v0, v1;
    unpack(elt.first, v0, v1);
    uint32_t e = elt.second - n_vtx;
    m.e2vtx[2*e] = v0;
    m.e2vtx[2*e+1] = v1;
  }
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
  #if USE_P2
  compute_transport_P2(this, T);
  #else
  compute_transport_P1(this, T);
  #endif
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
  for (uint32_t i=0; i<N; ++i)
    psi.data[i] *= -1;

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
