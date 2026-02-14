#include <stdint.h>

#include "array.h"
#include "mesh.h"

/* Edge to Triangle adjacency table */
struct ETAdjacency
{
	struct ETri
	{
		Edge next;
		Edge prev;
		uint32_t tri;
	};
	TArray<uint32_t> degree;
	TArray<uint32_t> offset;
	TArray<ETri> etri;
	ETAdjacency(const Mesh &m);

	/* Auxiliary functions */
	void set_degree(const Mesh &m);
	void set_offset(size_t ne);
	void set_etri(const Mesh &m);
};
