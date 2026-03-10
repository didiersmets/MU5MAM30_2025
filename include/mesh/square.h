#ifndef SQUARE_H
#define SQUARE_H

#include <stdint.h>
#include "mesh.h"

/* Build a square grid mesh with n subdivisions per side */
int load_square_grid(Mesh &m, uint32_t nx, uint32_t ny);

#endif