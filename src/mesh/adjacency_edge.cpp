#include "mesh.h"
#include "adjacency_edge.h"
#include <cstdint>
#include <iostream>
using namespace std;

ETAdjacency::ETAdjacency(const Mesh &m)
{
	size_t ne = m.edge_count();
	size_t nt = m.triangle_count();

	/* Reserve memory */
	degree.resize(ne);
	offset.resize(ne);
	etri.resize(3 * nt);

	/* Build degree */
	set_degree(m);

	/* Build offset */
	set_offset(ne);

	/* Build etri */
	set_etri(m);
}

/* Auxiliary functions */

void ETAdjacency::set_degree(const Mesh &m)
{
	degree.set(0);
	size_t nt = m.triangle_count();

	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = m.indices[3 * tri];
		uint32_t b = m.indices[3 * tri + 1];
		uint32_t c = m.indices[3 * tri + 2];

		Edge ab(a, b), bc(b, c), ca(c, a);

		degree[*m.edge_idx.get(ab)]++;
		degree[*m.edge_idx.get(bc)]++;
		degree[*m.edge_idx.get(ca)]++;
	}
}

void ETAdjacency::set_offset(size_t ne)
{
	offset[0] = 0;
	for (size_t k = 1; k < ne; k++)
		offset[k] = offset[k - 1] + degree[k - 1];
}

void ETAdjacency::set_etri(const Mesh &m)
{
	ETri zero;
	zero.next.v0 = UINT32_MAX;
	zero.prev.v0 = UINT32_MAX;
	etri.set(zero);

	size_t nt = m.triangle_count();

	ETri current_pair;
	size_t k;

	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t tri_verts[3] = {m.indices[3 * tri + 0], m.indices[3 * tri + 1], m.indices[3 * tri + 2]};

		for (size_t t = 0; t < 3; t++)
		{
			uint32_t a = tri_verts[t];
			uint32_t b = tri_verts[(t + 1) % 3];
			uint32_t c = tri_verts[(t + 2) % 3];

			Edge ab = Edge(a, b);
			Edge bc = Edge(b, c);
			Edge ca = Edge(c, a);

			uint32_t ab_idx = *m.edge_idx.get(ab);

			k = offset[ab_idx];
			while (k < offset[ab_idx] + degree[ab_idx] && etri[k].next.v0 != UINT32_MAX)
				k++;

			if (k < offset[ab_idx] + degree[ab_idx])
			{
				etri[k].next = bc;
				etri[k].prev = ca;
				etri[k].tri = tri;
			}
		}
	}
}
