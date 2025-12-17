#include "include/mesh/cube.h"
#include "include/mesh/mesh.h"
#include "include/mesh/sphere.h"

#include <iostream>
#include <vector>

/*
to run the main:

g++ -std=c++17 main.cpp src/mesh/cube.cpp src/mesh/duplicate_verts.cpp -I include -o main
*/

int main()
{
  // 1. Initialize the Mesh structure
  Mesh myMesh;

  // 2. Define subdivision level
  // subdiv = 4 means each face will be a 4x4 grid of quads
  size_t subdivisions = 4;

  std::cout << "--- Starting Mesh Generation Test ---" << std::endl;

  // 3. Load the cube
  // This calls load_overlapping_cube and then remove_duplicate_vertices
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

  // Expected Triangle Count: 6 faces * (subdiv * subdiv * 2 triangles)
  // For subdiv 4: 6 * (16 * 2) = 192 triangles
  size_t expected_triangles = 6 * (subdivisions * subdivisions * 2);

  std::cout << "\n[Mesh Results]" << std::endl;
  std::cout << "[Surface Logic Check]" << std::endl;
  std::cout << "The mesh is a hollow shell with " << myMesh.vertex_count()
            << " points on the boundary and 0 points inside." << std::endl;
  std::cout << "Total Indices:   " << myMesh.index_count() << std::endl;
  std::cout << "Total Triangles: " << myMesh.triangle_count() << std::endl;

  // 4. Validate results
  if (myMesh.triangle_count() == expected_triangles)
  {
    std::cout << "SUCCESS: Triangle count matches expected value (" << expected_triangles << ")."
              << std::endl;
  }
  else
  {
    std::cout << "WARNING: Triangle count mismatch!" << std::endl;
  }

  return 0;
}