#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

#include <iostream>
using namespace std;

static void load_cube_vertices(Vec3 *pos, size_t subdiv);
static void load_cube_indices(uint32_t *idx, size_t subdiv);

void rotate_y_to_x(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos)
{
	size_t nb_points = ref_pos.size;
	for (size_t i = 0; i < nb_points; i++)
	{
		pos[i].x = ref_pos[i].y;
		pos[i].y = ref_pos[i].x;
		pos[i].z = ref_pos[i].z;
	}
}

void rotate_y_to_z(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos)
{
	size_t nb_points = ref_pos.size;
	for (size_t i = 0; i < nb_points; i++)
	{
		pos[i].x = ref_pos[i].x;
		pos[i].y = ref_pos[i].z;
		pos[i].z = ref_pos[i].y;
	}
}

void flip_my_to_py(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos)
{
	size_t nb_points = ref_pos.size;
	for (size_t i = 0; i < nb_points; i++)
	{
		pos[i].x = ref_pos[i].x;
		pos[i].y = -ref_pos[i].y;
		pos[i].z = ref_pos[i].z;
	}
}

void flip_mx_to_px(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos)
{
	size_t nb_points = ref_pos.size;
	for (size_t i = 0; i < nb_points; i++)
	{
		pos[i].x = -ref_pos[i].x;
		pos[i].y = ref_pos[i].y;
		pos[i].z = ref_pos[i].z;
	}
}

void flip_mz_to_pz(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos)
{
	size_t nb_points = ref_pos.size;
	for (size_t i = 0; i < nb_points; i++)
	{
		pos[i].x = ref_pos[i].x;
		pos[i].y = ref_pos[i].y;
		pos[i].z = -ref_pos[i].z;
	}
}

int load_cube(Mesh &m, size_t subdiv)
{
	/* Check subdiv is reasonable and return error if not */
	if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */)
	{
		return (-1);
	}

	/* Assume subdiv corresponds to the number of cells in one direction */
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
	size_t nb_points = subdiv + 1;
	double radius = 1.;
	size_t nb_cells = subdiv;

	/* Build the positions of one reference face (arbitrary face y = -1) */
	TArray<Vec3> face_y_m1(pow(nb_points, 2));
	for (size_t i = 0; i < nb_points; i++)
		for (size_t j = 0; j < nb_points; j++)
		{
			face_y_m1[nb_points * i + j].x = radius * (-1. + 2. * j / nb_cells);
			face_y_m1[nb_points * i + j].y = -1. * radius;
			face_y_m1[nb_points * i + j].z = radius * (-1. + 2. * i / nb_cells);
		}

	/* Generate all the other faces by geometric transformations */
	TArray<Vec3> face_x_m1(pow(nb_points, 2));
	TArray<Vec3> face_z_m1(pow(nb_points, 2));
	TArray<Vec3> face_y_1(pow(nb_points, 2));
	TArray<Vec3> face_x_1(pow(nb_points, 2));
	TArray<Vec3> face_z_1(pow(nb_points, 2));

	rotate_y_to_x(face_y_m1, face_x_m1);
	rotate_y_to_z(face_y_m1, face_z_m1);
	flip_my_to_py(face_y_m1, face_y_1);
	flip_mx_to_px(face_x_m1, face_x_1);
	flip_mz_to_pz(face_z_m1, face_z_1);

	size_t index;

	for (size_t i = 0; i < 6; i++)
		for (size_t j = 0; j < pow(nb_points, 2); j++)
		{
			switch (i)
			{
			case 0:
				index = j;
				pos[index] = face_y_m1[j];
				break;
			case 1:
				index = size_t(pow(nb_points, 2)) + j;
				pos[index] = face_x_m1[j];
				break;
			case 2:
				index = size_t(pow(nb_points, 2)) * 2 + j;
				pos[index] = face_z_m1[j];
				break;
			case 3:
				index = size_t(pow(nb_points, 2)) * 3 + j;
				pos[index] = face_y_1[j];
				break;
			case 4:
				index = size_t(pow(nb_points, 2)) * 4 + j;
				pos[index] = face_x_1[j];
				break;
			default:
				index = size_t(pow(nb_points, 2)) * 5 + j;
				pos[index] = face_z_1[j];
				break;
			}
		}
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	size_t nb_cells = subdiv;
	size_t nb_points = subdiv + 1;

	/* Build the indices of one reference face */
	size_t current_index = 0;
	TArray<uint32_t> idx_y_m1(6 * pow(nb_cells, 2));
	for (size_t i = 0; i < nb_cells; i++)
		for (size_t j = 0; j < nb_cells; j++)
		{
			current_index = 6 * (nb_cells * i + j);

			/* Lower-left triangle of the cell */
			idx_y_m1[current_index] = nb_points * i + j;
			idx_y_m1[current_index + 1] = nb_points * i + j + 1;
			idx_y_m1[current_index + 2] = nb_points * (i + 1) + j;

			/* Upper-right triangle of the cell */
			idx_y_m1[current_index + 3] = nb_points * i + j + 1;
			idx_y_m1[current_index + 4] = nb_points * (i + 1) + j + 1;
			idx_y_m1[current_index + 5] = nb_points * (i + 1) + j;
		}

	/* Generate all the indices by translating the ones of the reference face */
	size_t nb_indices = 6 * pow(nb_cells, 2);
	for (size_t i = 0; i < 6; i++)
		for (size_t j = 0; j < nb_indices; j++)
			idx[i * nb_indices + j] = i * pow(nb_points, 2) + idx_y_m1[j];
}
