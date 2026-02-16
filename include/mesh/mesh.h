#pragma once
#include <stdint.h>

#include "array.h"
#include "vec3.h"

#include <map>
#include <utility>

struct Mesh {
    TArray<Vec3> positions;
    TArray<uint32_t> indices;
    TArray<float> attr;

    size_t vertex_count() const { return positions.size; }
    size_t index_count() const { return indices.size; }
    size_t triangle_count() const { return indices.size / 3; }
};

// DECLARACIÓN de la función (afuera del struct)
void build_edge_numbering(const Mesh &m,
                          std::map<std::pair<uint32_t,uint32_t>, uint32_t> &edge_id,
                          uint32_t &edge_count);
