#pragma once
#include "duplicate_verts.h"
#include "mesh.h"

int load_cube(Mesh& m, size_t subdiv);

/**
 * @brief Generates a cube mesh with overlapping vertices at edges/corners.
 *        Each face is generated separately, creating duplicate vertices.
 *        Triangles are generated with diagonal from top-left to bottom-right,
 *        all with normals pointing outward from the cube.
 *
 * @param m The mesh to fill
 * @param subdiv The number of triangles for each face of the cube is: 2×subdiv^2.
 * @return int vertices per face
 */

template <typename T>
int load_overlapping_cube(Mesh& m, size_t subdiv);

/**
 * @brief Generates a face of the cube.
 * @param invert: determines the winding order (triangle orientation).
 */
template <typename T, bool invert>
void build_cube_face(Mesh&    m,
                     size_t   subdiv,
                     TVec3<T> p_corner,
                     TVec3<T> dir_x,
                     TVec3<T> dir_y,
                     TVec3<T> normal)
{
  double h                = 2.0 / subdiv;
  size_t start_idx        = m.positions.size;
  size_t vertices_per_row = subdiv + 1;

  // 1. Generate vertices for this face
  for (size_t row = 0; row <= subdiv; row++)
  {
    for (size_t col = 0; col <= subdiv; col++)
    {
      TVec3<T> p = p_corner + dir_x * (T) (col * h) + dir_y * (T) (row * h);
      // Normalization of p can be added here if needed for the sphere
      m.positions.push_back(p);  // not normalized for the cube
    }
  }

  // 2. Generate triangles with specific winding
  for (size_t row = 0; row < subdiv; row++)
  {
    for (size_t col = 0; col < subdiv; col++)
    {
      uint32_t tl = (uint32_t) (start_idx + row * vertices_per_row + col);
      uint32_t tr = tl + 1;
      uint32_t bl = (uint32_t) (start_idx + (row + 1) * vertices_per_row + col);
      uint32_t br = bl + 1;

      if constexpr (invert)
      {
        m.indices.push_back(tl);
        m.indices.push_back(br);
        m.indices.push_back(tr);
        m.indices.push_back(tl);
        m.indices.push_back(bl);
        m.indices.push_back(br);
      }
      else
      {
        m.indices.push_back(tl);
        m.indices.push_back(tr);
        m.indices.push_back(br);
        m.indices.push_back(tl);
        m.indices.push_back(br);
        m.indices.push_back(bl);
      }
    }
  }
}

template <typename T>
int load_overlapping_cube(Mesh& m, size_t subdiv)
{
  m.positions.clear();
  m.indices.clear();

  // Call build_cube_face 6 times with appropriate corners and directions
  build_cube_face<T, false>(m,
                            subdiv,
                            TVec3<T>(-1, 1, -1),
                            TVec3<T>(1, 0, 0),
                            TVec3<T>(0, 0, 1),
                            TVec3<T>(0, 1, 0));  // TOP
  build_cube_face<T, false>(m,
                            subdiv,
                            TVec3<T>(-1, -1, 1),
                            TVec3<T>(1, 0, 0),
                            TVec3<T>(0, 0, -1),
                            TVec3<T>(0, -1, 0));  // BOTTOM
  build_cube_face<T, true>(m,
                           subdiv,
                           TVec3<T>(-1, -1, 1),
                           TVec3<T>(1, 0, 0),
                           TVec3<T>(0, 1, 0),
                           TVec3<T>(0, 0, 1));  // FRONT
  build_cube_face<T, true>(m,
                           subdiv,
                           TVec3<T>(1, -1, -1),
                           TVec3<T>(-1, 0, 0),
                           TVec3<T>(0, 1, 0),
                           TVec3<T>(0, 0, -1));  // BACK
  build_cube_face<T, false>(m,
                            subdiv,
                            TVec3<T>(1, -1, -1),
                            TVec3<T>(0, 0, 1),
                            TVec3<T>(0, 1, 0),
                            TVec3<T>(1, 0, 0));  // RIGHT
  build_cube_face<T, false>(m,
                            subdiv,
                            TVec3<T>(-1, -1, 1),
                            TVec3<T>(0, 0, -1),
                            TVec3<T>(0, 1, 0),
                            TVec3<T>(-1, 0, 0));  // LEFT

  return (int) m.vertex_count();
}