#include "boundary.h"
#include <unordered_map>
#include <algorithm>

static inline uint64_t edge_key(uint32_t a, uint32_t b)
{
    uint32_t lo = std::min(a,b);
    uint32_t hi = std::max(a,b);
    return (uint64_t(lo) << 32) | uint64_t(hi);//key=lo*2^32+hi, this way edge_key(a,b) and edge_key(b,a) are the same
}

void compute_boundary_vertices(const Mesh &m, TArray<uint8_t> &is_bnd)
{
    size_t N = m.vertex_count();
    size_t T = m.triangle_count();

    is_bnd.resize(N);
    for (size_t i = 0; i < N; ++i) is_bnd[i] = 0;
    std::unordered_map<uint64_t, uint32_t> edge_count;
    edge_count.reserve(3*T);

    const uint32_t *idx = m.indices.data;

    // count how many times each edge appears
    for (size_t t = 0; t < T; ++t)
    {
        uint32_t a = idx[3*t + 0];
        uint32_t b = idx[3*t + 1];
        uint32_t c = idx[3*t + 2];

        edge_count[edge_key(a,b)]++;
        edge_count[edge_key(b,c)]++;
        edge_count[edge_key(c,a)]++;
    }

    // Edge with count=1 are boundary edges, mark their vertices as boundary vertices
    for (auto &kv : edge_count)
    {
        if (kv.second == 1)
        {
           uint32_t a = uint32_t(kv.first >> 32);
            uint32_t b = uint32_t(kv.first & 0xffffffffu);
            is_bnd[a] = 1;
            is_bnd[b] = 1;
        }
    }
}
