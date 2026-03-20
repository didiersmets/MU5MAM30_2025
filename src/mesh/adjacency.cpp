#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m) 
	: degree(m.vertex_count(), 0), offset(m.vertex_count()), vtri(3*m.triangle_count())
{
	// Computing the degrees
	// Cost: O(n)
	for ( int i = 0; i < m.indices.size; i++ ) {
		this->degree[m.indices[i]]++;
	}
	// Computing the offsets
	// Cost: O(n)
	this->offset[0] = 0;
	for ( uint32_t index = 1; index < m.vertex_count(); index++ ) {
		this->offset[index] = this->offset[index-1] + this->degree[index-1];
	}

	// Tracks the first available position for storing next adiacent triangle to a certain vertex
	TArray<uint32_t> filling(m.vertex_count(), 0);

	// Computing the vtri
	// Cost: O(n)
	for ( int i = 0; i < m.index_count(); i += 3 ) {
		uint32_t indexA = m.indices[i];
		uint32_t indexB = m.indices[i+1];
		uint32_t indexC = m.indices[i+2];

		this->vtri[this->offset[indexA]+filling[indexA]] = {indexB, indexC};
		this->vtri[this->offset[indexB]+filling[indexB]] = {indexC, indexA};
		this->vtri[this->offset[indexC]+filling[indexC]] = {indexA, indexB};

		filling[indexA]++;
		filling[indexB]++;
		filling[indexC]++;
	}
}
