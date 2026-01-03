#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
    size_t num_vertices = m.vertex_count();
    size_t num_triangles = m.triangle_count();
    
    // Step 1: Initialize degree array (count triangles per vertex)
    degree.resize(num_vertices);
    
    // Count how many triangles each vertex is part of
    for (size_t t = 0; t < num_triangles; t++) {
        uint32_t v0 = m.indices[3*t + 0];
        uint32_t v1 = m.indices[3*t + 1];
        uint32_t v2 = m.indices[3*t + 2];
        
        degree[v0]++;
        degree[v1]++;
        degree[v2]++;
    }
    
    // Step 2: Build offset array (cumulative sum of degrees)
    offset.resize(num_vertices);
    offset[0] = 0;
    for (size_t i = 1; i < num_vertices; i++) {
        offset[i] = offset[i-1] + degree[i-1];
    }
    
    // Step 3: Allocate vtri array
    size_t total_entries = 3 * num_triangles;
    vtri.resize(total_entries);
    
    // Step 4: Fill vtri array
    // Track current write position for each vertex
    TArray<uint32_t> current_pos(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        current_pos[i] = offset[i];
    }
    
    // For each triangle, add the (next, prev) pairs to each vertex
    for (size_t t = 0; t < num_triangles; t++) {
        uint32_t v0 = m.indices[3*t + 0];
        uint32_t v1 = m.indices[3*t + 1];
        uint32_t v2 = m.indices[3*t + 2];
        
        // For vertex v0: next=v1, prev=v2
        vtri[current_pos[v0]].next = v1;
        vtri[current_pos[v0]].prev = v2;
        current_pos[v0]++;
        
        // For vertex v1: next=v2, prev=v0
        vtri[current_pos[v1]].next = v2;
        vtri[current_pos[v1]].prev = v0;
        current_pos[v1]++;
        
        // For vertex v2: next=v0, prev=v1
        vtri[current_pos[v2]].next = v0;
        vtri[current_pos[v2]].prev = v1;
        current_pos[v2]++;
    }
}