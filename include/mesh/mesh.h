#pragma once
#include <stdint.h>

#include "array.h"
#include "vec3.h"

struct Mesh {
	TArray<Vec3> positions; // per-vertex positions
	TArray<uint32_t> indices; // triangle indices, each group of 3 consecutive uint32_t represents a triangle by indexing into the positions array
	TArray<float> attr; // per-vertex scalar attribute, used for visualization, deformation, etc. The value of attr[i] is associated to the vertex at positions[i]
	TArray<uint32_t> boundary; // explicit list of boundary vertex indices
	size_t vertex_count() const { return positions.size; }
	size_t index_count() const { return indices.size; }
	size_t triangle_count() const { return indices.size / 3; }
	size_t boundary_count() const { return boundary.size; }
};
