#include "mesh.h"
#include "adjacency.h"
#include <vector>

VTAdjacency::VTAdjacency(const Mesh &m)
{
	/* Your implementation goes here */
    const size_t nv = m.vertex_count();
    const size_t nt = m.triangle_count();
    const size_t total_entries = 3 * nt;

    degree.resize(nv);
    offset.resize(nv);
    vtri.resize(total_entries);

    for (size_t i = 0; i < nv; ++i) {
        degree[i] = 0;
        offset[i] = 0;
    }

    for (size_t k = 0; k < total_entries; ++k) {
        degree[m.indices[k]]++;
    }

    for (size_t v = 1; v < nv; ++v) {
        offset[v] = offset[v - 1] + degree[v - 1];
    }

    std::vector<size_t> cursor(nv);
    for (size_t v = 0; v < nv; ++v) {
        cursor[v] = offset[v];
    }

    // Remplir vtri
    for (size_t t = 0; t < nt; ++t) {
        uint32_t a = m.indices[3 * t];
        uint32_t b = m.indices[3 * t + 1];
        uint32_t c = m.indices[3 * t + 2];

        vtri[cursor[a]++] = {b, c};
        vtri[cursor[b]++] = {c, a};
        vtri[cursor[c]++] = {a, b};
    }

    assert(offset[nv - 1] + degree[nv - 1] == total_entries);
}
