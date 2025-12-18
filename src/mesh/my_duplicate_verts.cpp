#include <stddef.h>
#include <stdint.h>

#include "mesh.h"
#include "vec3.h"

#include <iostream>
#include <algorithm>
#include <unordered_map>
using namespace std;


size_t build_position_remap(const TArray<Vec3> *pos, size_t count, uint32_t *remap)
{
    unordered_map<Vec3, uint32_t, Vec3Hash> last_seen;    
    size_t unique_vertices = 0;

    for (size_t i = 0; i < count; ++i)
    {
        auto it = last_seen.find((*pos)[i]);
        if (it != last_seen.end())
            remap[i] = last_seen[(*pos)[i]]; 
        else
            remap[i] = unique_vertices++;
        last_seen[(*pos)[i]] = remap[i];
    }

    return unique_vertices;
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
    for (size_t i = 0; i < old_vertex_count; i++)
        new_pos[remap[i]] = m.positions[i];

    /* Fill the new indices array */
    for (size_t j = 0; j < index_count; j++)
        new_idx[j] = remap[m.indices[j]];

    m.positions.resize(new_vertex_count);
    for (size_t i = 0; i < new_vertex_count; ++i)
        m.positions[i] = new_pos[i];

    for (size_t i = 0; i < index_count; ++i)
        m.indices[i] = new_idx[i];
}