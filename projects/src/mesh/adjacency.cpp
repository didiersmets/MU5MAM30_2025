#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	/* Your implementation goes here */
	size_t nb_vertex = m.vertex_count();
	size_t nb_triangle = m.triangle_count();
	
	degree.resize(nb_vertex);
	for (size_t i = 0; i < nb_vertex; i++){
		degree[i] = 0;
	}
	for (size_t i = 0; i < nb_triangle*3; i++) {
		degree[m.indices[i]]++;
	}
	
	offset.resize(nb_vertex);
	offset[0]  = 0;
	for (size_t i = 1; i < nb_vertex; i++){
		offset[i] = offset[i-1] + degree[i-1];
	}

	vtri.resize(3*nb_triangle);
	for (size_t i = 0; i < nb_triangle; i++){
		uint32_t a = m.indices[3*i];
		uint32_t b = m.indices[3*i + 1];
		uint32_t c = m.indices[3*i + 2];
		vtri[offset[a]++] = { b, c };
		vtri[offset[b]++] = { c, a };
		vtri[offset[c]++] = { a, b };
	}

	for (size_t v = 0; v < nb_vertex; v++) {
		offset[v] -= degree[v];
	}
}
