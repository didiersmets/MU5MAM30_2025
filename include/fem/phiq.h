#pragma once

#include "vec3.h"

struct phiq {
    double phi[3]; // Shape function values at the triangle vertices
    double q[3];   // Gradient of shape functions
};