#pragma once
#include <stdint.h>

#include "array.h"
#include "vec3.h"

struct Mesh {
	TArray<Vec3> positions;
	TArray<uint32_t> indices;/*it maps the vertices to the triangles (every 3 consecutive indexes I obtain a triangle)*/
	TArray<float> attr;
	size_t vertex_count() const { return positions.size; } /*number of vertices*/
	size_t index_count() const { return indices.size; }
	size_t triangle_count() const { return indices.size / 3; }
};
