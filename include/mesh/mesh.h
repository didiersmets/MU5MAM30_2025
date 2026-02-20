#pragma once
#include <stdint.h>
#include <string>

#include "array.h"
#include "vec3.h"

using namespace std;

struct Mesh {
	TArray<Vec3f> positions;
	TArray<uint32_t> indices;
	TArray<float> attr;
	size_t vertex_count() const { return positions.size; }
	size_t index_count() const { return indices.size; }
	size_t triangle_count() const { return indices.size / 3; }
};

void save_to_obj(Mesh &m, string file_name);
