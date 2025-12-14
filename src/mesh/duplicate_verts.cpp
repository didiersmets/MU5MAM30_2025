#include "duplicate_verts.h"
#include <memory>

size_t build_position_remap(const Vec3 *pos, size_t count, uint32_t *remap) {
    int count_remapped = 0;
    for ( size_t i = 0; i < count; i++ ) {
        for (size_t j = 0; j <= i; j++ ) {
            if ( pos[i] == pos[j] ) {
                remap[i] = j;  // better performance if I do pos[j] - count_remapped[j]?
                count_remapped++;
                break;
            }
        }
    }
    return count - count_remapped;
}

void remove_duplicate_vertices(Mesh &m) {
    uint32_t *remap = (uint32_t*) std::malloc(m.vertex_count() * sizeof(uint32_t));
    build_position_remap(m.positions.data, m.vertex_count(), remap);
    uint32_t vertex_count = m.vertex_count();
    uint32_t index_count  = m.index_count();



    int cur_count = 0;
    for ( uint32_t i = 0; i < vertex_count; i++ ) {
        if ( i == remap[i] ) {
            m.positions[cur_count] = m.positions[i];
            for ( uint32_t j = 0; j < index_count; j++ ) {
                if ( remap[m.indices[j]] == i ) {
                    m.indices[j] = cur_count;
                }
            }
            cur_count++;
        }
    }

    m.positions.resize(cur_count);

    free(remap);
}