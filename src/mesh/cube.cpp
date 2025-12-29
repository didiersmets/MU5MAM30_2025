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

	/*
	we first compute the offset between vertices in each face, having 6 faces
	0 -> front
	1 -> back
	2 -> left
	3 -> right
	4 -> bottom
	5 -> top
	*/ 
	size_t face_offset[6];
	for( int i = 0; i < 6; i++ )
		face_offset[i] = i * POW2(subdiv + 1); //each face contains (subdiv+1)^2 vertices
	
	/*
	the idea is the following, given an edge of the cube and a number of segments (subdiv)
	we will end up having (subdiv + 1) vertices on that edge and considering a evenly spaced
	mesh we want each vertex to be at a distance of 1.0 / subdiv from each other
	Since our cube goes from -1 to +1 the edge lenght is 2 so the dx = 2 / subdiv
	so we have that the position of vertex i on that edge is:
	xi = -1 + i * dx 
	*/
	float *x_i = (float *)malloc( (subdiv + 1) * sizeof(float) );
	float dx = 2.0f / subdiv;

	for( int i = 0; i <= subdiv; i++ )
		x_i[i] = -1.0f + i * dx;

	/*
	now we need to assign for each vertex it's position depending on the face
	basically each face is a matrix 2D matrix cause we fix one of the 3 dimensions
	front -> fix y = -1 -> vary x and z
	back  -> fix y = +1 -> vary x and z
	and so on for each face
	*/
	for( int x = 0; x <= subdiv; x++){
		for( int y = 0; y <= subdiv; y++){
			pos[ face_offset[0]++ ] = { x_i[x], -1.0f, x_i[subdiv - y] }; 			 //front
			pos[ face_offset[1]++ ] = { x_i[subdiv - x], 1.0f, x_i[y] };			 //back
			
			pos[ face_offset[2]++ ] = { -1.0f, x_i[x], x_i[y] }; 					//left
			pos[ face_offset[3]++ ] = { 1.0f, x_i[subdiv - x], x_i[subdiv - y] };   //right

			pos[ face_offset[4]++ ] = { x_i[x], x_i[subdiv - y], -1.0f };			//bottom
			pos[ face_offset[5]++ ] = { x_i[x], x_i[y], 1.0f }; 					//top
		}
	}
	
	free(x_i);
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */

	/*
	what we want to do here is to create the triangulation, meaning that given the grid of
	vertices we want to have as output the triangle indexes that tells which 3 vertixes create
	each triangle
	so we loop for each face and for each square in the grid we create 2 triangles by "dividing"
	the square along the diagonal and we then save the indexes of the triangles
	given 4 vertexes:
	tri1 = { base , base + 1,						base + (num_verts_per_row) }
	tri2 = { base , base + (num_verts_per_row) + 1, base + (num_verts_per_row) }

	num_verts_per_row = subdiv + 1
	*/
	uint32_t step = 0;
	for( int face = 0; face < 6; face++){
		int offset = face * POW2(subdiv + 1); 
		for( int i = 0; i < subdiv; i++){
			for( int j = 0; j < subdiv; j++){
				uint32_t firstvert = offset + i * (subdiv + 1) + j;
				//first triangle
				idx[step++] = firstvert;
				idx[step++] = firstvert + 1;
				idx[step++] = firstvert + (subdiv + 1);
				//second triangle
				idx[step++] = firstvert;
				idx[step++] = firstvert + (subdiv + 1) + 1;
				idx[step++] = firstvert + (subdiv + 1);
			}
		}
	}
}
