#pragma once
#include <stdint.h>

#include "common/array.h"
#include "common/vec3.h"

struct Mesh {
	TArray<Vec3> positions;
	TArray<uint32_t> indices;
	TArray<float> attr;
	size_t vertex_count() const { return positions.size; }
	size_t index_count() const { return indices.size; }
	size_t triangle_count() const { return indices.size / 3; }
};
