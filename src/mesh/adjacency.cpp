#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	size_t vtx_count = m.vertex_count();
	size_t idx_count = 3*m.triangle_count();
	degree.resize(vtx_count);
	offset.resize(vtx_count);
	vtri.resize(idx_count);


	for (uint32_t i = 0; i <vtx_count; i ++){
		degree[m.indices[i]] += 1;
	};
	for (uint32_t i = 1; i < vtx_count; i ++){
		offset[i] += (offset[i-1] + degree[i-1]);
	};

	// count how many of the incident triangles to a certain vertex we already saw
	TArray<uint32_t> tricount (vtx_count, 0); 

	for(uint32_t i = 0; i < idx_count; i +=3){
		for(int j = 0; j < 3; j ++){
			uint32_t a = m.indices[i + j];
            uint32_t b = m.indices[i + (j + 1) % 3];
            uint32_t c = m.indices[i + (j + 2) % 3];

            uint32_t pos = offset[a] + tricount[a];
            vtri[pos].next = b;
            vtri[pos].prev = c;
            tricount[a]++;	
		};
	};

}
