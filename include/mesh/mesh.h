#pragma once
#include "common/array.h"
#include "common/vec3.h"

#include <stdint.h>
struct Mesh
{
  // Array of unique vertex positions in 3D space (Vertex Buffer)
  TArray<Vec3> positions;

  // Array of vertex indices defining triangle connectivity (Index Buffer)
  // Every 3 consecutive indices represent one triangle
  TArray<uint32_t> indices;

  // Additional vertex data (UVs, Normals, Colors, etc.)
  TArray<float> attr;

  // Returns the total number of unique vertices
  size_t vertex_count() const { return positions.size; }

  // Returns the total number of entries in the index buffer
  size_t index_count() const { return indices.size; }

  // Returns the total number of triangles (assuming a Triangle List primitive)
  size_t triangle_count() const { return indices.size / 3; }
};