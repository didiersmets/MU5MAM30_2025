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
	
	uint32_t n = subdiv + 1;

	TArray<float> coord_pos(n);

	for(uint32_t i = 0; i < n; i ++){
		coord_pos[i] = -1.0f + 2.0f * i / subdiv;
	}

	
	for(uint32_t i = 0; i < n; i ++){
		for(uint32_t j = 0; j < n; j ++){

		pos[n*i + j ] = {1.0, coord_pos[i], coord_pos[j] };
		pos[POW2(n) + n*i + j ] ={-1.0, coord_pos[subdiv - i], coord_pos[j] };

		pos[2*POW2(n) + n*i + j ] = {coord_pos[subdiv - i], 1.0, coord_pos[j] };
		pos[3*POW2(n) + n*i + j ] = {coord_pos[i], -1.0, coord_pos[j]};

		pos[4*POW2(n) +n*i + j ] = {coord_pos[i], coord_pos[j], 1.0 };
		pos[5*POW2(n) +n*i + j ] = {coord_pos[subdiv-i], coord_pos[j], -1.0 };
		}
	}
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	
	uint32_t n = subdiv+1;

	uint32_t idx_counter = 0; // count how many triangles we have assigned
	
	for(int f = 0; f < 6; f ++) { //face we are currently looking at
		for(uint32_t i = 0; i < n-1; i ++){ // 'cloumn' on face
			for(uint32_t j = 0; j < n-1; j ++){ // 'row' on face 

			uint32_t v0, v1, v2, v3;
			v0 = f*POW2(n) + i + n*j;
			v1 = f*POW2(n) + i +1 +n*j;
			v2 = f*POW2(n) + i + 1+ n*(j +1);
			v3 = f*POW2(n) + i +  n*(j +1);

			idx[idx_counter] = v0;
			idx[idx_counter+1] = v1;
			idx[idx_counter+2] = v2;
			idx_counter += 3;

			idx[idx_counter] = v0;
			idx[idx_counter + 1] = v2;
			idx[idx_counter +2] = v3;

			idx_counter += 3;

			}
		}
	}
}
