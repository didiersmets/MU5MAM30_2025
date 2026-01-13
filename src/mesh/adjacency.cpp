#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
  /* (done) Your implementation goes here */

  size_t n_vtx = m.vertex_count();
  size_t size_ind = m.index_count();
  size_t n_tri = m.triangle_count();

  /* Create degree array */
  degree.resize(n_vtx);
  for (uint32_t v=0; v<n_vtx; ++v)
    degree[v] = 0;
  for (size_t i=0; i<size_ind; ++i) {
    uint32_t v = m.indices[i];
    ++degree[v];
  }

  /* Create offset array */
  offset.resize(n_vtx);
  offset[0] = 0;
  for (uint32_t k=1; k<n_vtx; ++k)
    offset[k] = offset[k-1] + degree[k-1];

  /* Create vtri array */
  vtri.resize(size_ind);
  TArray<uint32_t> loc_offset(n_vtx, 0);

  for (size_t tri=0; tri<n_tri; ++tri) {
    for (uint32_t a=0; a<3; ++a) {
      uint32_t b = ((a < 2) ? a+1 : 0);
      uint32_t c = ((a >= 1) ? a-1 : 2);

      uint32_t va = m.indices[3*tri+a];
      uint32_t vb = m.indices[3*tri+b];
      uint32_t vc = m.indices[3*tri+c];

      uint32_t insert_pos_a = offset[va] + loc_offset[va];
      vtri[insert_pos_a].next = vb;
      vtri[insert_pos_a].prev = vc;

      ++loc_offset[va];
    }
  }
}
