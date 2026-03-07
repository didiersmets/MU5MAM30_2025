#include "mesh.h"
#include "adjacency.h"
#include <cstdint>
#include <iostream>
using namespace std;

VTAdjacency::VTAdjacency(const Mesh &m)
{
	if (!m.is_periodic)
	{
		size_t nv = m.vertex_count();
		size_t nt = m.triangle_count();

		/* Reserve memory */
		degree.resize(nv);
		offset.resize(nv);
		vtri.resize(3 * nt);

		/* Build degree */
		set_degree(m.indices.data, nt);

		/* Build offset */
		set_offset(nv);

		/* Build vtri */
		set_vtri(m.indices.data, nt);
	}
	else
	{
		size_t nb_dofs = m.periodic_dofs_count;
		size_t nt = m.triangle_count();

		/* Reserve memory */
		degree.resize(nb_dofs);
		offset.resize(nb_dofs);
		vtri.resize(3 * nt);

		/* Build degree */
		set_degree_per(m);

		/* Build offset */
		set_offset_per(m);

		/* Build vtri */
		set_vtri_per(m);
	}
}

/* Auxiliary functions */

void VTAdjacency::set_degree(uint32_t *idx, size_t nt)
{
	degree.set(0);
	for (size_t tri = 0; tri < nt; tri++)
		for (size_t l = 0; l < 3; l++)
			degree[idx[3 * tri + l]] += 1;
}
void VTAdjacency::set_degree_per(const Mesh &m)
{
	degree.set(0);
	size_t nt = m.triangle_count();
	for (size_t tri = 0; tri < nt; tri++)
		for (size_t l = 0; l < 3; l++)
			degree[m.dof_map[m.indices.data[3 * tri + l]]] += 1;
}

void VTAdjacency::set_offset(size_t nv)
{
	offset[0] = 0;
	for (size_t k = 1; k < nv; k++)
		offset[k] = offset[k - 1] + degree[k - 1];
}

void VTAdjacency::set_offset_per(const Mesh &m)
{
	offset[0] = 0;
	size_t nb_dofs = m.periodic_dofs_count;
	for (size_t k = 1; k < nb_dofs; k++)
		offset[k] = offset[k - 1] + degree[k - 1];
}

void VTAdjacency::set_vtri(uint32_t *idx, size_t nt)
{
	VTri zero;
	zero.next = UINT32_MAX;
	zero.prev = UINT32_MAX;
	vtri.set(zero);

	VTri current_pair;
	size_t k;
	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t tri_verts[3] = {idx[3 * tri + 0], idx[3 * tri + 1], idx[3 * tri + 2]};

		for (size_t t = 0; t < 3; t++)
		{
			uint32_t a = tri_verts[t];
			uint32_t b = tri_verts[(t + 1) % 3];
			uint32_t c = tri_verts[(t + 2) % 3];

			k = offset[a];
			while (k < offset[a] + degree[a] && vtri[k].next != UINT32_MAX)
				k++;

			if (k < offset[a] + degree[a])
			{
				current_pair.next = b;
				current_pair.prev = c;
				vtri[k] = current_pair;
			}
		}
	}
}

void VTAdjacency::set_vtri_per(const Mesh &m)
{
	size_t nt = m.triangle_count();

	VTri zero;
	zero.next = UINT32_MAX;
	zero.prev = UINT32_MAX;
	vtri.set(zero);

	VTri current_pair;
	size_t k;
	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t tri_verts[3] = {m.dof_map[m.indices.data[3 * tri + 0]],
								 m.dof_map[m.indices.data[3 * tri + 1]],
								 m.dof_map[m.indices.data[3 * tri + 2]]};

		for (size_t t = 0; t < 3; t++)
		{
			uint32_t a = tri_verts[t];
			uint32_t b = tri_verts[(t + 1) % 3];
			uint32_t c = tri_verts[(t + 2) % 3];

			k = offset[a];
			while (k < offset[a] + degree[a] && vtri[k].next != UINT32_MAX)
				k++;

			if (k < offset[a] + degree[a])
			{
				current_pair.next = b;
				current_pair.prev = c;
				vtri[k] = current_pair;
			}
		}
	}
}
