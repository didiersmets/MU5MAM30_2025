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
	size_t pas { 2/subdiv };
	size_t n { subdiv +1 };


	int k { 0 };


	for (float a = -1; a <= 1 ; a+pas)
	{
		for (float b = -1; b <= 1; b+pas)
		{
			//face z=1 :
			Vec3 v {b, a, 1};
			pos[k] = v;

			//face y = 1
			v = {b, 1, a};
			pos[k + n*n] = v;

			//face x = 1
			v = {1.,b,a};
			pos[k + 2*n*n] = v;

			//face z = -1
			v = {b, a, -1.};
			pos[k + 3*n*n] = v;

			//face y = -1
			v = {b,-1.,a};
			pos[k + 4*n*n] = v;

			//face x = -1
			v = {-1., b, a};
			pos[k + 5*n*n] = v;

			k++;
		}
	}

	return;
}



static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */
	size_t n {subdiv+1};
	size_t tri_per_face {6*subdiv*subdiv};

	for (int i = 0; i < subdiv; ++i)
	{
		for (int j = 0; j < subdiv; ++j)
		{
			int k {6*(3*i + j)};

			//triangles on the face z=1
			idx[k] = n*i + j;
			idx[k+1] = n*i + j+1;
			idx[k+2] = n*(i+1) + j;

			idx[k+3] = n*i + j+1;
			idx[k+4] = n*(i+1) + j+1;
			idx[k+5] = n*(i+1) + j+1;


			//triangles on the face y=1
			idx[k+tri_per_face] = (n*i + j) + n*n;
			idx[k+1+tri_per_face] = n*i + j+1 + n*n;
			idx[k+2+tri_per_face] = n*(i+1) + j + n*n;

			idx[k+3+tri_per_face] = n*i + j+1 + n*n;
			idx[k+4+tri_per_face] = n*(i+1) + j+1 + n*n;
			idx[k+5+tri_per_face] = n*(i+1) + j+1 + n*n;


			//triangles on the face x=1
			idx[k+2*tri_per_face] = (n*i + j) + 2*n*n;
			idx[k+1+2*tri_per_face] = n*i + j+1 + 2*n*n;
			idx[k+2+2*tri_per_face] = n*(i+1) + j + 2*n*n;

			idx[k+3+2*tri_per_face] = n*i + j+1 + 2*n*n;
			idx[k+4+2*tri_per_face] = n*(i+1) + j+1 + 2*n*n;
			idx[k+5+2*tri_per_face] = n*(i+1) + j+1 + 2*n*n;


			//triangles on the face z=-1
			idx[k+3*tri_per_face] = (n*i + j) + 3*n*n;
			idx[k+1+3*tri_per_face] = n*i + j+1 + 3*n*n;
			idx[k+2+3*tri_per_face] = n*(i+1) + j + 3*n*n;

			idx[k+3+3*tri_per_face] = n*i + j+1 + 3*n*n;
			idx[k+4+3*tri_per_face] = n*(i+1) + j+1 + 3*n*n;
			idx[k+5+3*tri_per_face] = n*(i+1) + j+1 + 3*n*n;


			//triangles on the face y=-1
			idx[k+4*tri_per_face] = (n*i + j) + 4*n*n;
			idx[k+1+4*tri_per_face] = n*i + j+1 + 4*n*n;
			idx[k+2+4*tri_per_face] = n*(i+1) + j + 4*n*n;

			idx[k+3+4*tri_per_face] = n*i + j+1 + 4*n*n;
			idx[k+4+4*tri_per_face] = n*(i+1) + j+1 + 4*n*n;
			idx[k+5+4*tri_per_face] = n*(i+1) + j+1 + 4*n*n;


			//triangles on the face x=-1
			idx[k+5*tri_per_face] = (n*i + j) + 5*n*n;
			idx[k+1+5*tri_per_face] = n*i + j+1 + 5*n*n;
			idx[k+2+5*tri_per_face] = n*(i+1) + j + 5*n*n;

			idx[k+3+5*tri_per_face] = n*i + j+1 + 5*n*n;
			idx[k+4+5*tri_per_face] = n*(i+1) + j+1 + 5*n*n;
			idx[k+5+5*tri_per_face] = n*(i+1) + j+1 + 5*n*n;

		}
	}
	
	return;
}
