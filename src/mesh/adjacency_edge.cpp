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

	/* Build vtri */
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
	size_t nt = m.triangle_count();
	uint32_t k = 0;
	uint32_t j = 0;

	/* Create a cursor to track the triangles that share a same edge */
	TArray<uint32_t> offset_cursor(offset.size);
	for (size_t i = 0; i < offset.size; ++i)
		offset_cursor[i] = offset[i];

	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = m.indices[3 * tri];
		uint32_t b = m.indices[3 * tri + 1];
		uint32_t c = m.indices[3 * tri + 2];

		Edge ab = Edge(a, b);
		Edge bc = Edge(b, c);
		Edge ca = Edge(c, a);

		k = *m.edge_idx.get(ab);
		j = offset_cursor[k]++;
		etri[j].next = bc;
		etri[j].prev = ca;
		etri[j].tri = tri;

		/* Edge bc */
		k = *m.edge_idx.get(bc);
		j = offset_cursor[k]++;
		etri[j].next = ca;
		etri[j].prev = ab;
		etri[j].tri = tri;

		/* Edge ca */
		k = *m.edge_idx.get(ca);
		j = offset_cursor[k]++;
		etri[j].next = ab;
		etri[j].prev = bc;
		etri[j].tri = tri;
	}
}
