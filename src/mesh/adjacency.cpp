#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	size_t vtx_count = m.vertex_count();
	degree.resize(vtx_count);
	offset.resize(vtx_count);
	vtri.resize(3*m.triangle_count());

	for(size_t i = 0; i < vtx_count;++i){
		degree[i] = 0;
		offset[i] = 0;
	}
	//Count the degree of each vertex
	for(size_t i = 0; i < m.index_count(); ++i){
		degree[m.indices[i]]++;
	}
	//Compute the offset array
	offset[0] = 0;
	for(size_t i = 1; i < vtx_count; ++i){
		offset[i] = offset[i-1] + degree[i-1];
	}
	//Resize the vtri array
	vtri.resize(m.index_count());
	
	//Fill the vtri array
	for(size_t i = 0; i < m.triangle_count(); ++i){
		uint32_t a = m.indices[3*i];
		uint32_t b = m.indices[3*i + 1];
		uint32_t c = m.indices[3*i +2];
		//Vertex a
		vtri[offset[a]++].next = b;
		vtri[offset[a]++].prev = c;
		//Vertex b
		vtri[offset[b]++].next = c;
		vtri[offset[b]++].prev = a;
		//Vertex c
		vtri[offset[c]++].next = a;
		vtri[offset[c]++].prev = b;
	}

	//Restore the offset array
	for(size_t i = vtx_count - 1; i > 0; --i){
		offset[i] -= degree[i];
	}
}