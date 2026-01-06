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
	size_t n = subdiv + 1;
	float step = 2.0f / subdiv;
	
	// Generate 6 faces of the cube
	// Each face has n x n vertices
	size_t offset = 0;
	
	// Face 0: Front (z = +1)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			pos[offset++] = Vec3(-1.0f + j * step, -1.0f + i * step, 1.0f);
		}
	}
	
	// Face 1: Back (z = -1)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			pos[offset++] = Vec3(1.0f - j * step, -1.0f + i * step, -1.0f);
		}
	}
	
	// Face 2: Right (x = +1)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			pos[offset++] = Vec3(1.0f, -1.0f + i * step, -1.0f + j * step);
		}
	}
	
	// Face 3: Left (x = -1)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			pos[offset++] = Vec3(-1.0f, -1.0f + i * step, 1.0f - j * step);
		}
	}
	
	// Face 4: Top (y = +1)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			pos[offset++] = Vec3(-1.0f + j * step, 1.0f, 1.0f - i * step);
		}
	}
	
	// Face 5: Bottom (y = -1)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			pos[offset++] = Vec3(-1.0f + j * step, -1.0f, -1.0f + i * step);
		}
	}
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	size_t n = subdiv + 1;
	size_t idx_offset = 0;
	
	// Generate indices for all 6 faces
	for (size_t face = 0; face < 6; face++) {
		uint32_t face_offset = face * n * n;
		
		// For each quad in the subdivision grid
		for (size_t i = 0; i < subdiv; i++) {
			for (size_t j = 0; j < subdiv; j++) {
				// Vertex indices for the current quad
				uint32_t v0 = face_offset + i * n + j;
				uint32_t v1 = face_offset + i * n + (j + 1);
				uint32_t v2 = face_offset + (i + 1) * n + (j + 1);
				uint32_t v3 = face_offset + (i + 1) * n + j;
				
				// First triangle (v0, v1, v2)
				idx[idx_offset++] = v0;
				idx[idx_offset++] = v1;
				idx[idx_offset++] = v2;
				
				// Second triangle (v0, v2, v3)
				idx[idx_offset++] = v0;
				idx[idx_offset++] = v2;
				idx[idx_offset++] = v3;
			}
		}
	}
}
