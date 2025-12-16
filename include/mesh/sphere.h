#include "mesh.h"

// it follows the same strategy as the cube generation code, but when adding a new vertex, it
// normalizes its position vector to project it onto the sphere surface.
int load_sphere(Mesh& m, size_t subdiv);
