#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	/* Your implementation goes here */
	const size_t nv = m.vertex_count();
	const size_t nt = m.triangle_count();

	// initialize arrays
	degree.resize(nv);
	for (size_t i = 0; i < nv; i++)
		degree[i] = 0u;

	offset.resize(nv);
	for (size_t i = 0; i < nv; i++)
		offset[i] = 0u;

	vtri.resize(3 * nt);
	for (size_t i = 0; i < 3 * nt; i++) {
		vtri[i].next = 0u;
		vtri[i].prev = 0u;
	}

	// count degrees (number of incident triangles per vertex)
	for (size_t t = 0; t < nt; t++) {
		const uint32_t a = m.indices[3 * t + 0];
		const uint32_t b = m.indices[3 * t + 1];
		const uint32_t c = m.indices[3 * t + 2];
		degree[a]++;
		degree[b]++;
		degree[c]++;
	}

	// build offsets (prefix sums)
	uint32_t acc = 0u;
	for (size_t i = 0; i < nv; i++) {
		offset[i] = acc;
		acc += degree[i];
	}

	// temporary cursor to keep insertion positions per vertex
	TArray<uint32_t> cur(nv, 0u);

	// fill vtri: for triangle (a,b,c) we add entries for a:(next=b,prev=c),
	// b:(next=c,prev=a), c:(next=a,prev=b)
	for (size_t t = 0; t < nt; t++) {
		const uint32_t a = m.indices[3 * t + 0];
		const uint32_t b = m.indices[3 * t + 1];
		const uint32_t c = m.indices[3 * t + 2];

		uint32_t pos = offset[a] + cur[a]++;
		vtri[pos].next = b;
		vtri[pos].prev = c;

		pos = offset[b] + cur[b]++;
		vtri[pos].next = c;
		vtri[pos].prev = a;

		pos = offset[c] + cur[c]++;
		vtri[pos].next = a;
		vtri[pos].prev = b;
	}
}

