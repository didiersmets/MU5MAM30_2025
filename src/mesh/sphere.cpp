// same as cube.cpp but for sphere
#include "mesh/sphere.h"

#include "mesh/cube.h"
/*
int load_sphere(Mesh& m, size_t subdiv)
{
  load_overlapping_sphere<float>(m, subdiv);
  remove_duplicate_vertices(m);
  return 0;
}
*/

int load_sphere(Mesh& m, size_t subdiv)
{
  if (int res = load_cube(m, subdiv))
    return (res);

  Vec3*  pos       = m.positions.data;
  size_t vtx_count = m.positions.size;
  for (size_t i = 0; i < vtx_count; ++i)
  {
    pos[i] = normalized(pos[i]);
    // pos[i] *= 5.0f;  // scale sphere radius to 5.0
  }

  return (0);
}