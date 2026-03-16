#pragma once
#include <stdint.h>

#include "array.h"
#include "vec3.h"

struct Mesh {
	TArray<Vec3> positions;
	TArray<uint32_t> indices;
	TArray<float> attr;

	TArray<uint8_t> boundary; //新增:标记边界点的数组,0表示内部点;1表示边界点
	
	size_t vertex_count() const { return positions.size; }
	size_t index_count() const { return indices.size; }
	size_t triangle_count() const { return indices.size / 3; }
};
