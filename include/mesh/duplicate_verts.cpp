#include "duplicate_verts.h"
#include <cstdlib> 

size_t build_position_remap(const Vec3 *pos, size_t count, uint32_t *remap)
{
      if(count == 0){
        return 0;
      }

      size_t unique_count = 0;

      for(size_t i = 0; i < count; i++){
       // assume unique until proven otherwise
        remap[i] = (uint32_t)i;

        for(size_t j=0; j < i; j++){

            // duplicate vertices check
            if(pos[i] == pos[j]){

                // point to the same unique vertex
                  remap[i] = remap[j]; 
                  break;
               }
        }
        // no duplicates found
        if(remap[i] == i){
            unique_count ++;
        }
    }
        
    return unique_count;
      
}


void remove_duplicate_vertices(Mesh& m)
{
    if (m.vertex_count() == 0) return;
    
    // Build remap
    uint32_t* remap = (uint32_t*)malloc(m.vertex_count() * sizeof(uint32_t));
    size_t unique_count = build_position_remap(m.positions.data, 
                                                m.vertex_count(), 
                                                remap);
    
    // Compact positions and update remap
    uint32_t write_pos = 0;
    for (uint32_t i = 0; i < m.vertex_count(); i++) {
        if (i == remap[i]) {
            m.positions[write_pos] = m.positions[i];
            remap[i] = write_pos;
            write_pos++;
        } else {
            remap[i] = remap[remap[i]];
        }
    }
    
    // Update indices
    for (uint32_t j = 0; j < m.index_count(); j++) {
        m.indices[j] = remap[m.indices[j]];
    }
    
    
    m.positions.resize(unique_count);
    free(remap);
}