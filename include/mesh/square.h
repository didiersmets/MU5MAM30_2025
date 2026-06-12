#pragma once

#include "mesh.h"

// Build a square mesh directly into `Mesh` (positions in XY plane, z=0)
void build_square_mesh(Mesh *m, size_t N, double a = 1.0, double b = 1.0);
