#pragma once
#include <stdint.h>

#include "array.h"
#include "vec3.h"
#include <string>

struct Mesh {
	TArray<Vec3> positions; //coordinates of verts
	TArray<uint32_t> indices; //every 3 indicies define a triangle
	TArray<float> attr; //extra data probalbly
	size_t vertex_count() const { return positions.size; }
	size_t index_count() const { return indices.size; }
	size_t triangle_count() const { return indices.size / 3; }
};

void save_to_obj(Mesh &m, std::string file_name);