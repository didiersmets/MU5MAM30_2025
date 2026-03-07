#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "torus.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

#include <iostream>
using namespace std;

static void load_torus_vertices(Vec3 *pos, size_t subdiv);
static void load_torus_indices(uint32_t *idx, size_t subdiv);
static uint32_t load_torus_periodic_map(uint32_t *per_map, uint32_t *dof_map, size_t subdiv);

int load_torus(Mesh &m, size_t subdiv)
{
    /* Check subdiv is reasonable and return error if not */
    if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */)
    {
        return (-1);
    }

    m.is_periodic = true;

    /* Assume subdiv corresponds to the number of cells in one direction */
    size_t n = subdiv + 1;

    /* Reserve memory for vertices and indices */
    m.positions.resize(POW2(n));
    m.indices.resize(6 * POW2(subdiv));

    m.periodic_map.resize(POW2(n));
    m.dof_map.resize(POW2(n));

    /* First build vertices of the torus surface  */
    /* See below for implementation */
    load_torus_vertices(m.positions.data, subdiv);

    /* Build corresponding triangulation indices */
    /* See below for implementation */
    load_torus_indices(m.indices.data, subdiv);

    /* Build the periodic map */
    m.periodic_dofs_count = load_torus_periodic_map(m.periodic_map.data, m.dof_map.data, subdiv);

    /* Finally attach faces between themselves */
    /* Implementation in src/duplicate_verts.cpp */
    remove_duplicate_vertices(m);


    return (0);
}

static void load_torus_vertices(Vec3 *pos, size_t subdiv)
{
    size_t nb_points = subdiv + 1;
    double radius_x = 1.;
    double radius_y = 1.;
    size_t nb_cells = subdiv;

    /* Build the positions of one reference face (arbitrary face z = 0) */
    TArray<Vec3> face_z_0(pow(nb_points, 2));
    for (size_t i = 0; i < nb_points; i++)
        for (size_t j = 0; j < nb_points; j++)
        {
            pos[nb_points * i + j].x = radius_x * (-1. + 2. * j / nb_cells);
            pos[nb_points * i + j].y = radius_y * (-1. + 2. * i / nb_cells);
            pos[nb_points * i + j].z = 0.;
        }
}

static void load_torus_indices(uint32_t *idx, size_t subdiv)
{
    size_t nb_cells = subdiv;
    size_t nb_points = subdiv + 1;

    /* Build the indices */
    size_t current_index = 0;
    for (size_t i = 0; i < nb_cells; i++)
        for (size_t j = 0; j < nb_cells; j++)
        {
            current_index = 6 * (nb_cells * i + j);

            /* Lower-left triangle of the cell */
            idx[current_index] = nb_points * i + j;
            idx[current_index + 1] = nb_points * i + j + 1;
            idx[current_index + 2] = nb_points * (i + 1) + j;

            /* Upper-right triangle of the cell */
            idx[current_index + 3] = nb_points * i + j + 1;
            idx[current_index + 4] = nb_points * (i + 1) + j + 1;
            idx[current_index + 5] = nb_points * (i + 1) + j;
        }
}

static uint32_t load_torus_periodic_map(uint32_t *per_map, uint32_t *dof_map, size_t subdiv)
{
    size_t nb_points = subdiv + 1;

    /* Build the indices */
    for (size_t i = 0; i < nb_points; i++)
    {
        for (size_t j = 0; j < nb_points; j++)
        {
            size_t ii = (i == nb_points - 1) ? 0 : i;
            size_t jj = (j == nb_points - 1) ? 0 : j;

            per_map[i * nb_points + j] = ii * nb_points + jj;
        }
    }

    uint32_t counter = 0;
    for (size_t i = 0; i < POW2(nb_points); ++i)
    {
        if (per_map[i] == i)
        {
            dof_map[i] = counter++;
        }
    }

    for (size_t i = 0; i < POW2(nb_points); ++i)
    {
        dof_map[i] = dof_map[per_map[i]];
    }

    return counter;
}
