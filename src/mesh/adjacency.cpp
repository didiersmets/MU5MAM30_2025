#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	/* Your implementation goes here */
	  /*degree-->Count how many triangles touch each vertex*/
	
	const size_t V = m.vertex_count();
    const size_t T = m.triangle_count();
	degree.resize(V, 0);/*initialize degree to 0 for each vertex: first entry give me the size, second entry is the value to initialize with*/
	for (size_t t=0;t<T;t++){ /*for each triangle*/
		uint32_t a = m.indices[3*t+0]; /*first vertex of triangle t*/
		uint32_t b = m.indices[3*t+1]; /*second vertex of triangle t*/
		uint32_t c = m.indices[3*t+2]; /*third vertex of triangle t*/
		degree[a]++; /*increment degree of vertex a*/
		degree[b]++; /*increment degree of vertex b*/
		degree[c]++; /*increment degree of vertex c*/
	}
	/*offset-->For each vertex v, offset[v] = index in vtri where the data for vertex v begins
	 se faccio vtri[offset[v]]posso accedere all'elemento dell'array vtri che mi permette di vedere quali sono i triangoli contenenti v */
	offset.resize(V,0);
	size_t sum=0;
	for(v=0;v<V;v++){
		if (v==0){
			offset[0]=0; /*first vertex always starts at position 0*/
		}else{
			offset[v]= offset[v-1]+degree[v-1]; /*offset of vertex v is offset of previous vertex + degree of previous vertex*/
	}
	sum+=degree[v];/*to compute the number of entryes of vtri*/
	}
	vtri.resize(sum); /*total number of entries in vtri is 3 times the number of triangles*/
	TArray<uint32_t> cursor = offset;/*copio il mio offset in un altro array che mi serve per tenere traccia di dove sono*/
	/*vtri-->For each vertex, list of triangles incident to that vertex*/
	 /*total number of entries in vtri is 3 times the number of triangles*/
	for (size_t i=0;i<T;i++){ /*total number of entries in vtri is 3 times the number of triangles*/
		uint32_t a= m.indices[3*i+0]; /*first vertex of triangle i*/
		uint32_t b= m.indices[3*i+1]; /*second vertex of triangle i*/
		uint32_t c= m.indices[3*i+2]; /*third vertex of triangle i*/
		/*For vertex a*/
		vtri[ cursor[a] ]= {0,0}; 
		cursor[a]++;/*since when I will have again vertex a my cursor moved on to the next cell containing the information of the next triangle*/
		/*For vertex b*/
		vtri[ cursor[b] ]= {0,0};	
		cursor[b]++;
		/*For vertex c*/
		vtri[ cursor[c] ]= {0,0};		
		cursor[c]++;
	}
	TArray<uint32_t> cursor = offset;
	for (size_t t=0;t<T;t++){
		uint32_t a = m.indices[3*t+0]; /*first vertex of triangle t*/
		uint32_t b = m.indices[3*t+1]; /*second vertex of triangle t*/
		uint32_t c = m.indices[3*t+2]; /*third vertex of triangle t*/
		/*For vertex a*/
		vtri[ cursor[a] ].next = b;
		cursor[a]++;/*since when I will have again vertex a my cursor moved on to the next cell containing the information of the next triangle*/
		/*For vertex b*/
		vtri[ cursor[b] ].next = c;
		vtri[ cursor[b] ].prev = a;
		cursor[b]++;
		/*For vertex c*/
		vtri[ cursor[c] ].next = a;
		vtri[ cursor[c] ].prev = b;
		cursor[c]++;
	}
	

}
