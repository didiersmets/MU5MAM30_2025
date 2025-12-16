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

#include "duplicate_verts.h"

#include <cstdlib>  // Required for malloc/free

size_t build_position_remap(const Vec3* pos, size_t count, uint32_t* remap)
{
  size_t unique_count = 0;

  for (size_t i = 0; i < count; i++)
  {
    bool found = false;
    // Search if this position has already been encountered
    for (size_t j = 0; j < i; j++)
    {
      if (pos[i] == pos[j])
      {
        remap[i] = remap[j];  // Map to the index of the first occurrence
        found    = true;
        break;
      }
    }

    // If not found, it's a unique vertex
    if (!found)
    {
      remap[i] = (uint32_t) i;
      unique_count++;
    }
  }
  return unique_count;
}

void remove_duplicate_vertices(Mesh& m)
{
  if (m.vertex_count() == 0)
    return;

  // Allocate memory for the remap table
  uint32_t* remap = (uint32_t*) std::malloc(m.vertex_count() * sizeof(uint32_t));

  // Step 1: Identify duplicates
  size_t unique_count = build_position_remap(m.positions.data, m.vertex_count(), remap);

  // Step 2: Create a mapping for the new condensed position array
  // We need to know where each old unique vertex will end up in the new array
  uint32_t* final_map      = (uint32_t*) std::malloc(m.vertex_count() * sizeof(uint32_t));
  uint32_t  cur_unique_idx = 0;

  for (uint32_t i = 0; i < m.vertex_count(); i++)
  {
    if (i == remap[i])
    {  // It's a first occurrence
      m.positions[cur_unique_idx] = m.positions[i];
      final_map[i]                = cur_unique_idx;
      cur_unique_idx++;
    }
    else
    {
      // It's a duplicate, point it to the already assigned new index
      final_map[i] = final_map[remap[i]];
    }
  }

  // Step 3: Update the index buffer to point to new locations
  for (uint32_t j = 0; j < m.index_count(); j++)
  {
    m.indices[j] = final_map[m.indices[j]];
  }

  // Step 4: Finalize
  m.positions.resize(unique_count);

  std::free(remap);
  std::free(final_map);
}
