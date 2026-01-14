#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	size_t nv = m.vertex_count();
	size_t index = m.index_count();
	size_t nt = m.triangle_count();
	int degree[nv];
	int offset[nv];
	VTri vtri[3*nt];
	
	for(int i=0;i<nv;i++){
		degree[i]=0;
		offset[i]=0;
	}
	for(int i=0;i<index;i++){
		degree[m.indices[index]]+=1;
	}
	for(int i=1;i<nv,i++){
		offset[i]=offset[i-1]+degree[i-1];
	}
	for(int i=0;i<nt,i++){
		int a = m.indices[3*i];
		int b = m.indices[3*i+1];
		int c = m.indices[3*i+2];
		vtri[offset[a]++]={b,c};
		vtri[offset[b]++]={c,a};
		vtri[offset[c]++]={a,b};
	}

	/* Your implementation goes here */
}
