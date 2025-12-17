#include "mesh.h"

/*
This structure stores only the "skin" of the object.
 * There are no internal vertices (volume is empty).
 * Ideal for surface-only PDEs (Laplace-Beltrami).
 */

//  Idead taken from cube.h , here the cube is projected onto a sphere.
// it follows the same strategy as the cube generation code,
// but when adding a new vertex,
// it normalizes its position vector to project it onto the sphere surface.
int load_sphere(Mesh& m, size_t subdiv);
