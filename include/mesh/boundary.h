#pragma once
#include "mesh.h"
#include "array.h"
#include <stdint.h>

// Restituisce un array is_bnd di dimensione N:
// is_bnd[i] = 1 se il vertice i è sul bordo
void compute_boundary_vertices(const Mesh &m, TArray<uint8_t> &is_bnd);