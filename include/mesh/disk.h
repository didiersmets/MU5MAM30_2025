#pragma once

#include "mesh.h"

// Build a disk mesh directly into `Mesh` (positions in XY plane, z=0)
void build_disk_mesh(Mesh *m, size_t N, double R = 1.0);
