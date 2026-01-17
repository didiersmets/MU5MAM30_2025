/*
build_position_remap: Identifies duplicate vertices and generates a mapping to unify them (Vertex
Welding).

remove_duplicate_vertices: Optimizes the mesh by removing redundant vertex data and updating indices
to point to unique locations.

workflow:

input: "dirty mesh" with duplicate vertices

positions: [A, B, C, C, B, D]
indices: [0, 1, 2, 3, 4, 5]

remap: [0, 1, 2, 2, 1, 3]
positions after removal: [A, B, C, D]
indices after update: [0, 1, 2, 2, 1, 3]
*/

#include "mesh/duplicate_verts.h"

#include "common/array.h"
#include "common/hash.h"
#include "common/hash_table.h"
#include "common/vec3.h"
#include "mesh/mesh.h"

#include <assert.h>
#include <cmath>  // Required for std::abs
#include <stddef.h>
#include <stdint.h>

// Helper to check if two vertices are close enough to be merged
static bool are_vertices_close(const Vec3& a, const Vec3& b)
{
  // Using a tolerance allows merging vertices that are mathematically
  // the same but different due to floating point rounding (e.g. 0.99999 vs 1.0)
  const double epsilon = 1e-5;
  return (std::abs(a.x - b.x) < epsilon) && (std::abs(a.y - b.y) < epsilon) &&
         (std::abs(a.z - b.z) < epsilon);
}

// O(N^2) implementation: Robust for Vertex Welding with Epsilon
size_t build_position_remap(const Vec3* pos, size_t count, uint32_t* remap)
{
  size_t unique_count = 0;

  for (size_t i = 0; i < count; ++i)
  {
    bool found = false;

    // Look backwards to see if this vertex already exists
    for (size_t j = 0; j < i; ++j)
    {
      if (are_vertices_close(pos[i], pos[j]))
      {
        remap[i] = remap[j];  // Reuse the index of the existing vertex
        found    = true;
        break;
      }
    }

    // If not found, it is a new unique vertex
    if (!found)
    {
      remap[i] = (uint32_t) i;  // Map to itself (initially)
      unique_count++;
    }
  }

  // Note: This logic assumes we will pack the array later.
  // The 'remap' now holds the index of the "first occurrence" for every vertex.
  // We need to compact this in the next step (in remove_duplicate_vertices).
  return unique_count;
}

void remove_duplicate_vertices(Mesh& m)
{
  // 1. Build the Remap Table
  TArray<uint32_t> remap(m.vertex_count());

  // Note: We perform the cast to const Vec3* here to match the signature
  build_position_remap(m.positions.data, m.vertex_count(), remap.data);

  // 2. Create the Condensed Position Array
  // We need to map the "old indices" to "new packed indices"
  TArray<uint32_t> old_to_new_map(m.vertex_count());
  size_t           new_idx_counter = 0;

  for (size_t i = 0; i < m.vertex_count(); ++i)
  {
    // If this vertex points to itself, it is a "Primary" vertex (the first of its kind)
    if (remap[i] == i)
    {
      // Move it to the new packed position
      m.positions[new_idx_counter] = m.positions[i];

      // Record where it went
      old_to_new_map[i] = (uint32_t) new_idx_counter;

      new_idx_counter++;
    }
  }

  // Fill in the map for the duplicate vertices
  for (size_t i = 0; i < m.vertex_count(); ++i)
  {
    if (remap[i] != i)
    {
      // If I am a duplicate of vertex J, my new index is the same as J's new index
      old_to_new_map[i] = old_to_new_map[remap[i]];
    }
  }

  // 3. Resize the vertex array to the correct size
  m.positions.resize(new_idx_counter);

  // 4. Update all triangles to point to the new indices
  for (size_t k = 0; k < m.indices.size; ++k)
  {
    m.indices[k] = old_to_new_map[m.indices[k]];
  }
}