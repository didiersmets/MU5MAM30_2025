#include "mesh.h"
#include "adjacency.h"
#include <cstdint>
#include <iostream>
using namespace std;

VTAdjacency::VTAdjacency(const Mesh &m)
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

/* Auxiliary functions */

void VTAdjacency::set_degree(uint32_t *idx, size_t nt)
{
	degree.set(0);
	for (size_t tri = 0; tri < nt; tri++)
		for (size_t l = 0; l < 3; l++)
			degree[idx[3 * tri + l]] += 1;
}

void VTAdjacency::set_offset(size_t nv)
{
	offset[0] = 0;
	for (size_t k = 1; k < nv; k++)
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
		/* Triplet (a, b, c) */
		k = offset[idx[3 * tri]];
		while (k < offset[idx[3 * tri]] + degree[idx[3 * tri]] && vtri[k].next != UINT32_MAX)
			k++;

		if (k < offset[idx[3 * tri]] + degree[idx[3 * tri]])
		{
			current_pair.next = idx[3 * tri + 1];
			current_pair.prev = idx[3 * tri + 2];
			vtri[k] = current_pair;
		}

		/* Triplet (b, c, a) */
		k = offset[idx[3 * tri + 1]];
		while (k < offset[idx[3 * tri + 1]] + degree[idx[3 * tri + 1]] && vtri[k].next != UINT32_MAX)
			k++;

		if (k < offset[idx[3 * tri + 1]] + degree[idx[3 * tri + 1]])
		{
			current_pair.next = idx[3 * tri + 2];
			current_pair.prev = idx[3 * tri];
			vtri[k] = current_pair;
		}

		/* Triplet (c, a, b) */
		k = offset[idx[3 * tri + 2]];
		while (k < offset[idx[3 * tri + 2]] + degree[idx[3 * tri + 2]] && vtri[k].next != UINT32_MAX)
			k++;

		if (k < offset[idx[3 * tri + 2]] + degree[idx[3 * tri + 2]])
		{
			current_pair.next = idx[3 * tri];
			current_pair.prev = idx[3 * tri + 1];
			vtri[k] = current_pair;
		}
	}
}
