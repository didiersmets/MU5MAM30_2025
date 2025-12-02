#include <stddef.h>
#include <stdint.h>

#include "mesh.h"
#include "vec3.h"

#include <iostream>
#include <algorithm>
using namespace std;

size_t build_position_remap(const TArray<Vec3> *pos, size_t count, uint32_t *remap)
{
    size_t j = -1;
    size_t erased = 0;
    for (size_t i = 0; i < count; i++)
    {
        remap[i] = i;
        j = 0;
        while (j < i && !((*pos)[i] == (*pos)[j]))
            j += 1;
        if (j != i)
        {
            remap[i] = j;
            erased += 1;
        }
    }

    return count - erased;
}

void remove_duplicate_vertices(Mesh &m)
{
    size_t old_vertex_count = m.vertex_count();
    size_t index_count = m.index_count();

    uint32_t *remap = new uint32_t[old_vertex_count];

    /* Build the positions remap */
    size_t new_vertex_count = build_position_remap(&m.positions, old_vertex_count, remap);

    TArray<Vec3> new_pos(new_vertex_count);
    TArray<uint32_t> new_idx(index_count);

    /* Fill the new positions array */
    size_t k = 0;
    for (size_t i = 0; i < old_vertex_count; i++)
        if (remap[i] == i)
        {
            new_pos[k] = m.positions[i];
            k += 1;
        }

    /* Correct the indices to take the new positions into account */
    for (size_t i = 0; i < old_vertex_count; i++)
    {
        size_t j = 0;
        while (!(m.positions[i] == new_pos[j]))
            j += 1;
        remap[i] = j;
    }

    /* Fill the new indices array */
    for (size_t j = 0; j < index_count; j++)
        new_idx[j] = remap[m.indices[j]];

    m.positions.resize(new_vertex_count);
    for (size_t i = 0; i < new_vertex_count; ++i)
        m.positions[i] = new_pos[i];

    for (size_t i = 0; i < index_count; ++i)
        m.indices[i] = new_idx[i];
}