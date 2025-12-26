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
	size_t face_offset[6];
	for(size_t f = 0; f < 6; ++f){
		face_offset[f] = f * POW2(n);
	}
	Vec3 face_directions[6] = { Vec3::XAxis, -Vec3::XAxis,
		                        Vec3::YAxis, -Vec3::YAxis,
								Vec3::ZAxis, -Vec3::ZAxis};
	for(size_t x = 0; x < n; ++x){
		for(size_t y = 0; y < n; ++y){
			pos[face_offset[0]++] = Vec3(1, (float)x / subdiv * 2 - 1, (float)y / subdiv * 2 - 1); //Front face
			pos[face_offset[1]++] = Vec3(-1, (float)(subdiv - x)/ subdiv * 2 - 1, (float)y / subdiv * 2 - 1); //Back face
			pos[face_offset[2]++] = Vec3((float)(subdiv - y) / subdiv * 2 - 1, -1, (float)x / subdiv * 2 - 1); //Left face
			pos[face_offset[3]++] = Vec3((float)y / subdiv * 2 - 1, 1, (float)(subdiv-x) / subdiv * 2 - 1); //Right face
			pos[face_offset[4]++] = Vec3((float)x / subdiv * 2 - 1, (float)(subdiv - y) / subdiv * 2 - 1, -1); //Bottom face
			pos[face_offset[5]++] = Vec3((float)(subdiv - x) / subdiv * 2 - 1, (float)y / subdiv * 2 - 1, 1); //Top face
		}
	}       
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */
	size_t n = subdiv + 1;
	//Build triangulation indices
	for(int f = 0; f < 6; ++f){
		size_t offset = f * POW2(n);
		for(size_t i = 0; i < subdiv; ++i){
			for(size_t j = 0; j < subdiv; ++j){
				uint32_t base = (uint32_t)(offset + i * n + j);
				*idx++ = base;
				*idx++ = base + 1;
				*idx++ = base + n;
				*idx++ = base + 1;
				*idx++ = base + 1 + n;
				*idx++ = base + 1;
			}
		}
	}
}
