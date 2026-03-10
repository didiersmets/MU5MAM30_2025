
#include "math_utils.h"
#include "mesh.h"
#include "half_sphere.h"
#include "vec3.h"
#include <vector>
#include <stdint.h>

int load_half_sphere(Mesh &m, size_t subdiv)
{
	if (int res = load_half_cube(m, subdiv))
		return res;

	Vec3 *pos = m.positions.data;
	size_t vtx_count = m.positions.size;
	for (size_t i = 0; i < vtx_count; ++i) {
		pos[i] = normalized(pos[i]);
	}

	return 0;
}