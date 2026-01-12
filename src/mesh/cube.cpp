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

	/* Computing usefull values*/
	size_t n = subdiv+1;
	float step = 1/((float) subdiv);
	Vec3 dx = {step, 0.0, 0.0};
	Vec3 dy = {0.0, step, 0.0};
	Vec3 dz = {0.0, 0.0, step};

	/* Face order don't matter, we choose a dice ordering */
	/* Front -> Top -> Right -> Left -> Bottom -> Back */

	/* For each face we choose the "down-left" corner (from a rotation of the front face) as the first vertex of the face*/
    Vec3 face_vec_offset;
	size_t face_index_offset;
	/* Front */
	face_vec_offset = {-0.5, -0.5, -0.5};
	face_index_offset = 0;
	for (size_t i=0; i < n; i++) {
		for (size_t j=0; j<n; j++) {
			pos[face_index_offset + i + j*n] = face_vec_offset + ((float) i)*dx + ((float) j)*dz; 
		}
	}
	
	/* Top */
	face_vec_offset = {-0.5, -0.5, 0.5};
	face_index_offset += POW2(n);
	for (size_t i=0; i < n; i++) {
		for (size_t j=0; j<n; j++) {
			pos[face_index_offset + i + j*n] = face_vec_offset + ((float) i)*dx + ((float) j)*dy; 
		}
	}

	/* Right */
	face_vec_offset = {0.5, -0.5, -0.5};
	face_index_offset += POW2(n);
	for (size_t i=0; i < n; i++) {
		for (size_t j=0; j<n; j++) {
			pos[face_index_offset + i + j*n] = face_vec_offset + ((float) i)*dy + ((float) j)*dz; 
		}
	}

	/* Left */
	face_vec_offset = {-0.5, 0.5, -0.5};
	face_index_offset += POW2(n);
	for (size_t i=0; i < n; i++) {
		for (size_t j=0; j<n; j++) {
			pos[face_index_offset + i + j*n] = face_vec_offset - ((float) i)*dy + ((float) j)*dz; 
		}
	}

	/* Bottom */
	face_vec_offset = {-0.5, 0.5, -0.5};
	face_index_offset += POW2(n);
	for (size_t i=0; i < n; i++) {
		for (size_t j=0; j<n; j++) {
			pos[face_index_offset + i + j*n] = face_vec_offset + ((float) i)*dx - ((float) j)*dy; 
		}
	}

	/* Back */
	face_vec_offset = {0.5, 0.5, -0.5};
	face_index_offset += POW2(n);
	for (size_t i=0; i < n; i++) {
		for (size_t j=0; j<n; j++) {
			pos[face_index_offset + i + j*n] = face_vec_offset - ((float) i)*dx + ((float) j)*dz; 
		}
	}
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */

	size_t n = subdiv+1;
	/* We build the array face by face following the same face ordering as before */
	/* Front -> Top -> Right -> Left -> Bottom -> Back */
	size_t face_index_offset = 0;
	size_t vertex_idx_ptr = 0; 

	for (size_t face_counter = 0; face_counter < 6; face_counter++) {
		for (size_t i = 0; i < subdiv; i++) {
			for (size_t j = 0; j < subdiv; j++) {
				/* Two triangles per base vertex  = 6 index per base vertex*/
				u_int32_t base_vertex_index = face_index_offset + i + j*n;
				idx[vertex_idx_ptr++] = base_vertex_index;
				idx[vertex_idx_ptr++] = base_vertex_index + 1;
				idx[vertex_idx_ptr++] = base_vertex_index + n;
				idx[vertex_idx_ptr++] = base_vertex_index + 1;
				idx[vertex_idx_ptr++] = base_vertex_index + n + 1;
				idx[vertex_idx_ptr++] = base_vertex_index + n;
			}
		}
		face_index_offset += POW2(n);
	}
}
