#include "common/array.h"
#include "matrix/sparse_matrix.h"
#include "mesh/mesh.h"

struct NavierStokesSolver
{
  NavierStokesSolver(const Mesh& m);
  const Mesh& m;
  size_t      N;    // DoF
  double      vol;  // Surface(m), used for insuring zero mean to omega and psi

  TArray<double> omega;
  TArray<double> Momega;
  TArray<double> psi;

  CSRPattern P;  // Pattern arrays
  CSRMatrix  S;  // Stiffness matrix
  CSRMatrix  M;  // Mass matrix

  TArray<double> r;   // current residue r = Mf - Su
  TArray<double> p;   // internal for cg
  TArray<double> Ap;  // internal for cg

  bool inited;  // Initialization computes first residue and error

  size_t iter_max = 500;
  double tol      = 1e-6;

  double t;

  void   set_zero_mean(double* V);
  size_t compute_stream_function();
  void   compute_transport(double* T);
  void   time_step_coriolis(double dt, double nu, double omega_earth);
  void   time_step(double dt, double nu);
  void   compute_transport_coriolis(double* T, double omega_earth);
};