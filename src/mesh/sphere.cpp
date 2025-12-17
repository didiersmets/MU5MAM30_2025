// same as cube.cpp but for sphere
#include "mesh/sphere.h"
int load_sphere(Mesh& m, size_t subdiv)
{
  // Uncomment when it will be implemeneted
  load_overlapping_sphere<float>(m, subdiv);
  remove_duplicate_vertices(m);
  return 1;
}