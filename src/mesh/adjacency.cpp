#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m) 
{

	int nv { m.vertex_count() };
	int nt { m.triangle_count() };

	for (int i = 0; i < nv; ++i)
	{
		degree[i] = 0;
		offset[i] = 0;
	}


	//first : we complete degree
	for (int i = 0; i < nt; ++i)
	{
		uint32_t a = m.indices[3*i];
		uint32_t b = m.indices[3*i+1];
		uint32_t c = m.indices[3*i+2];

		degree[a] += 1;
		degree[b] += 1;
		degree[c] += 1;
	}



	//Second : we complete offset
	for (int i = 1; i < nv; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			offset[i] += degree[j];
		}
	}


	/*Then : for vtri  
			 we create a list of shape 2x3xnt 
			 containing the previous and the next ertex 
			 of each vertex in each triangle
	*/

	TArray<uint32_t> pre_vtri (2*3*nt, 0);
	TArray<uint32_t> ind_for_vtri;

	for (int i = 0; i < nv; ++i)
	{
		ind_for_vtri[i] = 2*offset[i];
	}

	for (int i = 0; i < nt; ++i)
	{
		uint32_t a { m.indices[3*i] };
		uint32_t b { m.indices[3*i+1] };
		uint32_t c { m.indices[3*i+2] };

		pre_vtri[ind_for_vtri[a]] = b;
		pre_vtri[ind_for_vtri[a]+1] = c;
		ind_for_vtri[a] += 2;

		pre_vtri[ind_for_vtri[b]] = c;
		pre_vtri[ind_for_vtri[b]+1] = a;
		ind_for_vtri[b] += 2;

		pre_vtri[ind_for_vtri[c]] = a;
		pre_vtri[ind_for_vtri[c]+1] = b;
		ind_for_vtri[c] += 2;

	}
	

	for (int i = 0; i < 3*nt; ++i)
	{
		VTri v_i;
		v_i.next = pre_vtri[2*i];
		v_i.prev = pre_vtri[2*i+1];

		vtri[i] = v_i;
	}


	return;
}
