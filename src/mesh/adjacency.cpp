#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	/* Your implementation goes here */
	size_t vcount = m.vertex_count();
	size_t trcount = m.triangle_count();
	size_t total_vtris = trcount * 3;

	//check that the mesh indexes for triangles are valid
	assert( m.index_count() == total_vtris );

	//resize the arrays and initialize to 0
	degree.resize(vcount);
	offset.resize(vcount);
	vtri.resize(total_vtris);
	for (size_t i = 0; i < vcount; i++) {
		degree[i] = 0;
		offset[i] = 0;
	}

	/*
	compute the degree -> how many triangles each vertex is part of:
	m.indexes contains the vertex indexes for each triangle, so 
	for each time we meet an index we increment it's degree
	*/
	for( int i = 0; i < total_vtris; i++)
		degree[ m.indices[i] ]++;

	/*
	to help with the CSR structure we need to compute the offsets:
	offset[v] is the starting index in vtri[] for vertex v’s adjacency list
	*/
	for( int v = 1; v < vcount; v++ )
		offset[v] = offset[v-1] + degree[v-1];

	/*
	what we have left is to fill the vtri[] array:
	For each vertex, its adjacent triangle's vertexes
	are stored as a double linked list in an array
	we first copy offsets in a tmp array to avoid overwriting them
	the reason is that we need to increment the positions as we add entries
	to the vtri array
	*/
	TArray<uint32_t> tmp_offset;
	tmp_offset.resize(vcount);
	for( int v = 0; v < vcount; v++)
		tmp_offset[v] = offset[v];

	for( int tr = 0; tr < trcount; tr++){
		//we get the 3 vertex indexes for the triangle tr
		uint32_t v0 = m.indices[tr * 3     ];
		uint32_t v1 = m.indices[tr * 3 + 1 ];
		uint32_t v2 = m.indices[tr * 3 + 2 ];

		//for each vertex we put the other two in its vtri struct
		vtri[ tmp_offset[v0] ] = { v1, v2 };
		vtri[ tmp_offset[v1] ] = { v2, v0 };
		vtri[ tmp_offset[v2] ] = { v0, v1 };

		//increment the tmp_offset for the next triangle for each vertex
		tmp_offset[v0]++;
		tmp_offset[v1]++;
		tmp_offset[v2]++;
	}
}
