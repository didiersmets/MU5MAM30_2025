#include "include/fem/P1.h"
#include "include/mesh/cube.h"
#include "include/mesh/mesh.h"
#include "include/mesh/sphere.h"

#include <iostream>
#include <vector>

/*
to run the main:
g++ -std=c++17 main.cpp src/matrix/sparse_matrix.cpp src/mesh/adjacency.cpp src/fem/P1.cpp
src/mesh/sphere.cpp src/mesh/duplicate_verts.cpp -I include -o main
*/

int main()
{
  Mesh myMesh;

  // Define subdivision level
  // subdiv = 4 means each face will be a 4x4 grid of quads
  size_t subdivisions = 4;

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
  std::cout << "CSR Pattern built with " << pattern.cols << " non-zero entries (lower triangle)."
            << std::endl;
  std::cout << "[Step 2] Assembling Mass Matrix..." << std::endl;
  CSRMatrix massMatrix;
  build_P1_mass_matrix(myMesh, pattern, massMatrix);
  std::cout << "Mass Matrix assembled successfully." << std::endl;
  std::cout << "Matrix Rows:     " << massMatrix.rows << std::endl;
  std::cout << "Matrix NNZ:      " << massMatrix.nnz << std::endl;

  return 0;
}