#include "half_sphere.h"
#include "sphere.h"
#include "mesh.h"
#include "vec3.h"
#include <stdint.h>
#include <vector>

int load_half_sphere(Mesh &m, size_t subdiv)
{
    // load a full sphere
    Mesh full;
    if (load_sphere(full, subdiv))
        return -1;

    size_t full_vtx = full.positions.size;
    size_t full_tri = full.triangle_count();

    // only keeps the points that z >= 0 
    std::vector<uint32_t> old_to_new(full_vtx, UINT32_MAX);
    uint32_t new_idx = 0;
    for (size_t i = 0; i < full_vtx; ++i) {
        if (full.positions[i][2] >= 0.0f) {
            old_to_new[i] = new_idx++;
        }
    }

    std::vector<uint32_t> new_indices;
    for (size_t t = 0; t < full_tri; ++t) {
        uint32_t a = full.indices[3*t+0];
        uint32_t b = full.indices[3*t+1];
        uint32_t c = full.indices[3*t+2];
        if (old_to_new[a] != UINT32_MAX &&
            old_to_new[b] != UINT32_MAX &&
            old_to_new[c] != UINT32_MAX) {
            new_indices.push_back(old_to_new[a]);
            new_indices.push_back(old_to_new[b]);
            new_indices.push_back(old_to_new[c]);
        }
    }

    m.positions.resize(new_idx);
    for (size_t i = 0; i < full_vtx; ++i) {
        if (old_to_new[i] != UINT32_MAX)
            m.positions[old_to_new[i]] = full.positions[i];
    }

    m.indices.resize(new_indices.size());
    for (size_t i = 0; i < new_indices.size(); ++i)
        m.indices[i] = new_indices[i];

    return 0;
}