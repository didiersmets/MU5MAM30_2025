#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "torus.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

static void load_torus_vertices(Vec3 *pos, size_t subdiv);
static void load_torus_indices(uint32_t *idx, size_t subdiv);

int load_torus(Mesh &m, size_t subdiv)
{
	/* Check subdiv is reasonable and return error if not */
	if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */) {
		return (-1);
	}

	/* Reserve memory for vertices and indices */
	m.positions.resize(POW2(subdiv));
	m.indices.resize(6*POW2(subdiv));

	/* First build vertices */
	/* See below for implementation */
	load_torus_vertices(m.positions.data, subdiv);

	/* Build corresponding triangulation indices */
	load_torus_indices(m.indices.data, subdiv);


	return (0);
}

void load_torus_vertices(Vec3 *pos, size_t subdiv)
{
    float step = 1/((float) subdiv);
    Vec3 dx = {step, 0.0, 0.0};
    Vec3 dy = {0.0, step, 0.0};

    for (size_t i = 0; i < subdiv; i++) {
        for (size_t j = 0; j < subdiv; j++) {
            pos[i + j*subdiv] = ((float) i) * dx + ((float) j) * dy;
        } 
    }
}

void load_torus_indices(uint32_t *idx, size_t subdiv)
{
    size_t vertex_idx_ptr = 0;

    for (size_t i = 0; i < subdiv; i++) {
        for (size_t j = 0; j < subdiv; j++) {
				uint32_t next_i = (i+1) % subdiv;
                uint32_t next_j = (j+1) % subdiv;
				idx[vertex_idx_ptr++] = i + j * subdiv;
				idx[vertex_idx_ptr++] = next_i + j * subdiv;
				idx[vertex_idx_ptr++] = i + next_j * subdiv;
				idx[vertex_idx_ptr++] = next_i + j * subdiv;
				idx[vertex_idx_ptr++] = next_i + next_j * subdiv;
				idx[vertex_idx_ptr++] = i + next_j * subdiv;
        }
    }
}