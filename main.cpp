#include "include/fem/P1.h"
#include "include/linalg/conjugate_gradient.h"
#include "include/mesh/cube.h"
#include "include/mesh/mesh.h"
#include "include/mesh/sphere.h"

#include <iostream>
#include <vector>

/*
TO RUN THE MAIN use myMakefile
*/

int main()
{
  Mesh myMesh;

  // Define subdivision level
  // subdiv = 4 means each face will be a 4x4 grid of quads
  // size_t subdivisions = 4;

  size_t subdivisions = 10;

  std::cout << "--- Starting Mesh Generation Test ---" << std::endl;

  // Load the sphere mesh with specified subdivisions
  // This calls load_overlapping_sphere and then remove_duplicate_vertices
  if (load_sphere(myMesh, subdivisions))
  {
    std::cout << "Sphere successfully generated and optimized." << std::endl;
  }
  else
  {
    std::cerr << "Error: Failed to generate sphere mesh." << std::endl;
    return -1;
  }

  // --- Verification Logic ---
  size_t expected_triangles = 6 * (subdivisions * subdivisions * 2);

  std::cout << "\n[Mesh Results]" << std::endl;
  std::cout << "[Surface Logic Check]" << std::endl;
  std::cout << "The mesh is a hollow shell with " << myMesh.vertex_count()
            << " vertex on the boundary and 0 points inside." << std::endl;
  std::cout << "Total Indices:   " << myMesh.index_count() << std::endl;
  std::cout << "Total Triangles: " << myMesh.triangle_count() << std::endl;
  size_t total_vtris = myMesh.triangle_count() * 3;  // number of vertex-triangle connections
  std::cout << "Total connections:     " << total_vtris << std::endl;

  // Validate results
  if (myMesh.triangle_count() == expected_triangles)
  {
    std::cout << "SUCCESS: Triangle count matches expected value (" << expected_triangles << ")."
              << std::endl;
  }
  else
  {
    std::cout << "WARNING: Triangle count mismatch!" << std::endl;
  }

  CSRPattern pattern;
  build_P1_CSRPattern(myMesh, pattern);
  std::cout << "CSR Pattern built with " << std::endl;

  CSRMatrix massMatrix;
  CSRMatrix stiffMatrix;
  build_P1_mass_matrix(myMesh, pattern, massMatrix);
  build_P1_stiffness_matrix(myMesh, pattern, stiffMatrix);
  std::cout << "Mass Matrix (sym) assembled successfully." << std::endl;
  std::cout << "Matrix Rows:     " << massMatrix.rows << std::endl;
  std::cout << "Matrix NNZ:      " << massMatrix.nnz << std::endl;
  std::cout << "Stiffness Matrix (sym) assembled successfully." << std::endl;
  std::cout << "Matrix Rows:     " << stiffMatrix.rows << std::endl;
  std::cout << "Matrix NNZ:      " << stiffMatrix.nnz << std::endl;

  // Now massMatrix and stiffMatrix can be used for further FEM computations
  // A = M + S (SPD matrix)
  // b = M*f
  // solve A*x = b using conjugate gradient solver where x is the unknown DOFs

  CSRMatrix A;
  A.rows      = massMatrix.rows;
  A.cols      = massMatrix.cols;
  A.nnz       = massMatrix.nnz;
  A.symmetric = massMatrix.symmetric;
  A.row_start = massMatrix.row_start;
  A.col       = massMatrix.col;
  A.data.resize(A.nnz);

  // compute A = M + S
  for (size_t k = 0; k < A.nnz; k++)
  {
    A.data[k] = massMatrix.data[k] + stiffMatrix.data[k];
  }

  // --- SET UP TEST CASE
  // We choose u_exact = x^2 - z^2.
  // Since this is a spherical harmonic of degree k=2,
  // the source term f is given by: f = (-Delta + 1)u = 6u + u = 7u = 7(x^2 - z^2).

  size_t              n = A.cols;
  std::vector<double> f(n, 0.0);
  std::vector<double> u_exact(n, 0.0);

  for (size_t i = 0; i < n; ++i)
  {
    // Access vertex coordinates (adjust according to your Mesh class)
    double vx = myMesh.positions[i].x;
    // not used in this set up
    // double vy = myMesh.positions[i].y;
    double vz = myMesh.positions[i].z;

    // Define the exact solution and the source term
    u_exact[i] = vx * vx - vz * vz;
    f[i]       = 7.0 * (vx * vx - vz * vz);
  }

  // --- ASSEMBLE RIGHT-HAND SIDE (b = M * f) ---
  std::vector<double> b(n, 0.0);
  massMatrix.mvp(f.data(), b.data());

  // --- SOLVE THE SYSTEM ---
  std::vector<double> x(n, 0.0);  // Initial solution x0
  std::vector<double> r(n);
  std::vector<double> p(n);
  std::vector<double> Ap(n);

  double rel_error      = 0.0;
  double tolerance      = 1e-6;
  int    max_iterations = 1000;

  size_t iterations = conjugate_gradient_solve(A,
                                               b.data(),
                                               x.data(),
                                               r.data(),
                                               p.data(),
                                               Ap.data(),
                                               &rel_error,
                                               tolerance,
                                               max_iterations,
                                               false);

  std::cout << "CG converged in " << iterations << " iterations."
            << "the relative error is: " << rel_error << std::endl;
  return 0;
}