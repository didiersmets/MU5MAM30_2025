#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

static void load_cube_vertices(Vec3 *pos, size_t subdiv);
static void load_cube_indices(uint32_t *idx, size_t subdiv);

int load_cube(Mesh &m, size_t subdiv)
{
	/* Check subdiv is reasonable and return error if not */
	if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */) {
		return (-1);
	}

	size_t n = subdiv + 1;

	/* Reserve memory for vertices and indices */
	m.positions.resize(6 * POW2(n));
	m.indices.resize(36 * POW2(subdiv));

	/* First build vertices as six unattached faces of n^2 vertices each */
	/* See below for implementation */
	load_cube_vertices(m.positions.data, subdiv);

	/* Build corresponding triangulation indices */
	/* See below for implementation */
	load_cube_indices(m.indices.data, subdiv);

	/* Finally attach faces between themselves */
	/* Implementation in src/duplicate_verts.cpp */
	remove_duplicate_vertices(m);

	return (0);
}

static void load_cube_vertices(Vec3 *pos, size_t subdiv)
{
	/* Your implementation goes here */
	size_t n = subdiv + 1;
	size_t voff[6];
	for (int f = 0; f < 6; f++) {
		voff[f] = f * POW2(n);
	}
	/* Build grid direction */
	float *dir = static_cast<float *>(malloc(n * sizeof(float)));
	float *rev = static_cast<float *>(malloc(n * sizeof(float)));
	for (size_t i = 0; i < n; i++) {
		dir[i] = (2.0f * static_cast<float>(i) - static_cast<float>(subdiv)) / static_cast<float>(subdiv);
	}
	for (size_t i = 0; i < n; i++) {
		rev[i] = dir[subdiv - i];
	}
	for (size_t y = 0; y < n; y++) {
		for (size_t x = 0; x < n; x++) {
			pos[voff[0]++] = { dir[x], -1.0f, dir[y] }; /* Front  */
			pos[voff[1]++] = { rev[x], +1.0f, dir[y] }; /* Back   */
			pos[voff[2]++] = { -1.0f, rev[x], dir[y] }; /* Left   */
			pos[voff[3]++] = { +1.0f, dir[x], dir[y] }; /* Right  */
			pos[voff[4]++] = { dir[x], rev[y], -1.0f }; /* Bottom */
			pos[voff[5]++] = { dir[x], dir[y], +1.0f }; /* Top    */
		}
	}
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */
	size_t n = subdiv + 1;
	/* Build corresponding triangulation indices */
	for (int f = 0; f < 6; f++) {
		size_t offset = f * POW2(n);
		for (size_t i = 0; i < subdiv; i++) {
			for (size_t j = 0; j < subdiv; j++) {
				uint32_t base = static_cast<uint32_t>(i * n + j + offset);
				/* First tri */
				*idx++ = base;
				*idx++ = base + 1;
				*idx++ = base + 1 + n;
				/* Second tri */
				*idx++ = base;
				*idx++ = base + 1 + n;
				*idx++ = base + n;
			}
		}
	}
}

