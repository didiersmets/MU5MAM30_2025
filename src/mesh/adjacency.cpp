#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	/* Your implementation goes here */

	/* First we initialize degree, offset and vtri to the right size */
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	degree.resize(vtx_count);
	offset.resize(vtx_count);
	vtri.resize(3*tri_count);

	/* We compute degree[] by going through every indices and 
	adding one to degree[a] each time a is found in the indices*/

	uint32_t *idx = m.indices.data;
	size_t index_count = m.index_count();

	for (size_t i=0; i<index_count; i++) { 
		degree[idx[i]]++;
	}

	/* We compute offset once degree is fully filled */
	uint32_t running_sum = 0;
	for (size_t i=0; i<index_count; i++) {
		offset[i] = running_sum;
		running_sum += degree[i];
	}

	/* Finally we can compute vtri */
	/* We add a way to know how many incident triangles were found for a vertex during the loop*/
	TArray<uint32_t> known_tri;
	known_tri.resize(vtx_count);

	for (size_t index=0; index<index_count; index+=3) { /* we skip to the next triangle */
		uint32_t a_idx = idx[index];
		uint32_t b_idx = idx[index + 1];
		uint32_t c_idx = idx[index + 2];
		
		/* There are 3 incident triangle from this triangle, one for each of its vertex */
		VTri vtri_a = {b_idx, c_idx};
		VTri vtri_b = {c_idx, a_idx};
		VTri vtri_c = {a_idx, b_idx};

		/* We fill the vtri at the right spot in the array */
		vtri[offset[a_idx] + known_tri[a_idx]] = vtri_a;
		vtri[offset[b_idx] + known_tri[b_idx]] = vtri_b;
		vtri[offset[c_idx] + known_tri[c_idx]] = vtri_c;

		/* We update the known incident triangle counter for each vertex */
		known_tri[a_idx]++;
		known_tri[b_idx]++;
		known_tri[c_idx]++;
	}
}
