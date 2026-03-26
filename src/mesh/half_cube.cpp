#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "half_cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"


#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

static void load_half_cube_vertices(Vec3 *pos, size_t subdiv)
{
	size_t n = subdiv + 1;
	size_t k = 0;

	for (size_t j = 0; j < n; ++j) {
		float v = (float)j / (float)subdiv;          // [0,1]
		float z = v;                                 // [0,1]
		float y = -1.0f + 2.0f * v;                  // [-1,1]

		for (size_t i = 0; i < n; ++i) {
			float u = (float)i / (float)subdiv;
			float x = -1.0f + 2.0f * u;              // [-1,1]

			// 0) top: z = +1
			pos[k++] = Vec3{ x, y, 1.0f };
		}
	}

	for (size_t j = 0; j < n; ++j) {
		float v = (float)j / (float)subdiv;
		float z = v;                                 // [0,1]

		for (size_t i = 0; i < n; ++i) {
			float u = (float)i / (float)subdiv;
			float x = -1.0f + 2.0f * u;              // [-1,1]

			// 1) front: y = +1
			pos[k++] = Vec3{ x, 1.0f, z };
		}
	}

	for (size_t j = 0; j < n; ++j) {
		float v = (float)j / (float)subdiv;
		float z = v;                                 // [0,1]

		for (size_t i = 0; i < n; ++i) {
			float u = (float)i / (float)subdiv;
			float x = 1.0f - 2.0f * u;               // invertito per orientazione

			// 2) back: y = -1
			pos[k++] = Vec3{ x, -1.0f, z };
		}
	}

	for (size_t j = 0; j < n; ++j) {
		float v = (float)j / (float)subdiv;
		float z = v;                                 // [0,1]

		for (size_t i = 0; i < n; ++i) {
			float u = (float)i / (float)subdiv;
			float y = 1.0f - 2.0f * u;               // [+1,-1]

			// 3) right: x = +1
			pos[k++] = Vec3{ 1.0f, y, z };
		}
	}

	for (size_t j = 0; j < n; ++j) {
		float v = (float)j / (float)subdiv;
		float z = v;                                 // [0,1]

		for (size_t i = 0; i < n; ++i) {
			float u = (float)i / (float)subdiv;
			float y = -1.0f + 2.0f * u;              // [-1,1]

			// 4) left: x = -1
			pos[k++] = Vec3{ -1.0f, y, z };
		}
	}
}
static inline uint32_t local_vid(uint32_t i, uint32_t j, uint32_t n)
{
	return j * n + i;
}

static void add_face_indices(uint32_t *idx, uint32_t &k, uint32_t base, uint32_t subdiv)
{
	uint32_t n = subdiv + 1;

	for (uint32_t j = 0; j < subdiv; ++j) {
		for (uint32_t i = 0; i < subdiv; ++i) {
			uint32_t a = base + local_vid(i,   j,   n);
			uint32_t b = base + local_vid(i+1, j,   n);
			uint32_t c = base + local_vid(i+1, j+1, n);
			uint32_t d = base + local_vid(i,   j+1, n);

			idx[k++] = a;
			idx[k++] = b;
			idx[k++] = c;

			idx[k++] = a;
			idx[k++] = c;
			idx[k++] = d;
		}
	}
}

static void load_half_cube_indices(uint32_t *idx, size_t subdiv)
{
	uint32_t n = (uint32_t)subdiv + 1;
	uint32_t face_size = n * n;
	uint32_t k = 0;

	for (uint32_t f = 0; f < 5; ++f) {
		add_face_indices(idx, k, f * face_size, (uint32_t)subdiv);
	}
}

int load_half_cube(Mesh &m, size_t subdiv)
{
	if (subdiv <= 0 || subdiv > (1 << 14)) {
		return -1;
	}

	size_t n = subdiv + 1;

	// 5 facce invece di 6
	m.positions.resize(5 * POW2(n));
	m.indices.resize(5 * 6 * POW2(subdiv));

	load_half_cube_vertices(m.positions.data, subdiv);
	load_half_cube_indices(m.indices.data, subdiv);

	remove_duplicate_vertices(m);
	return 0;
}