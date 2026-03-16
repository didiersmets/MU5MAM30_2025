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
	const size_t n = subdiv + 1;

    size_t face_offset[6];
    for (size_t f = 0; f < 6; ++f) {
        face_offset[f] = f * POW2(n);
    }

    float *coord = (float *)safe_malloc(n * sizeof(float));
    for (size_t i = 0; i < n; ++i) {
        coord[i] = (2.0f * (float)i - (float)subdiv) / (float)subdiv;
    }

    // Remplissage des 6 faces
    for (size_t y = 0; y < n; ++y) {
        for (size_t x = 0; x < n; ++x) {
            const float dx = coord[x];
            const float dy = coord[y];
            const float rx = coord[subdiv - x];
            const float ry = coord[subdiv - y];

            pos[face_offset[0]++] = { dx, -1.0f, dy }; // Front
            pos[face_offset[1]++] = { rx,  1.0f, dy }; // Back
            pos[face_offset[2]++] = {-1.0f, rx,  dy }; // Left
            pos[face_offset[3]++] = { 1.0f, dx,  dy }; // Right
            pos[face_offset[4]++] = { dx,  ry, -1.0f}; // Bottom
            pos[face_offset[5]++] = { dx,  dy,  1.0f}; // Top
        }
    }

    free(coord);
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */
    const size_t n = subdiv + 1;
    size_t k = 0;

    for (size_t f = 0; f < 6; ++f) {
        const uint32_t face_offset = static_cast<uint32_t>(f * POW2(n));

        for (size_t row = 0; row < subdiv; ++row) {
            for (size_t col = 0; col < subdiv; ++col) {
                const uint32_t v0 = face_offset + static_cast<uint32_t>(row * n + col);
                const uint32_t v1 = v0 + 1;
                const uint32_t v2 = v0 + static_cast<uint32_t>(n);
                const uint32_t v3 = v2 + 1;

                // Triangle 1
                idx[k++] = v0;
                idx[k++] = v1;
                idx[k++] = v3;

                // Triangle 2
                idx[k++] = v0;
                idx[k++] = v3;
                idx[k++] = v2;
            }
        }
    }
}
