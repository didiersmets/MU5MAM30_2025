#pragma once

#include <stdint.h>

#include "array.h"
#include "my_mesh.h"

struct MyVTAdjacency {
	struct VTri {
		uint32_t next;
		uint32_t prev;
	};
	TArray<uint32_t> degree;
	TArray<uint32_t> offset;
	TArray<VTri> vtri;
	MyVTAdjacency(const MyMesh &m);
};
