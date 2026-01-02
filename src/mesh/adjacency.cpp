#include "mesh/adjacency.h"

#include "mesh/mesh.h"

VTAdjacency::VTAdjacency(const Mesh& m)
{
  /* Your implementation goes here */
  size_t vcount  = m.vertex_count();  // unique vertices in the mesh
  size_t trcount = m.triangle_count();
  // total number of indices in the mesh == (3 per triangle) * (number of triangles)
  size_t total_vtris = trcount * 3;

  // check that the mesh indexes for triangles are valid
  assert(m.index_count() == total_vtris);

  // resize the arrays and initialize to 0
  degree.resize(vcount);
  offset.resize(vcount);
  vtri.resize(total_vtris);
  for (int i = 0; i < vcount; i++)
  {
    degree[i] = 0;
    offset[i] = 0;
  }

  /*
  compute the degree for each vertex :
  counting how many triangles are connected to the specific vertex.
  m.indices[] contains all the vertex indexes for all triangles
  */
  for (int i = 0; i < total_vtris; i++)
    degree[m.indices[i]]++;

  /*
  to help with the CSR structure we need to compute the offsets:
  offset[v] is the starting index in vtri[] for vertex v’s adjacency list
  */
  for (int v = 1; v < vcount; v++)
    offset[v] = offset[v - 1] + degree[v - 1];

  TArray<uint32_t> tmp_offset;
  tmp_offset.resize(vcount);
  for (int v = 0; v < vcount; v++)
    tmp_offset[v] = offset[v];

  for (int tr = 0; tr < trcount; tr++)
  {
    // we get the 3 vertex indexes for the triangle tr
    uint32_t v0 = m.indices[tr * 3];
    uint32_t v1 = m.indices[tr * 3 + 1];
    uint32_t v2 = m.indices[tr * 3 + 2];

    // for each vertex of the triangle we add the other two vertices to its adjacency list
    vtri[tmp_offset[v0]] = {v1, v2};
    vtri[tmp_offset[v1]] = {v2, v0};
    vtri[tmp_offset[v2]] = {v0, v1};

    // increment the tmp_offset for the next triangle for each vertex
    tmp_offset[v0]++;
    tmp_offset[v1]++;
    tmp_offset[v2]++;
  }
  /* SUMMARY OF THE FINAL STRUCTURE (vtri):
   *
   * 1. BLOCK-BASED ORGANIZATION (CSR-STYLE)
   * The vtri array is a flat array (dimension: 3 * m.triangle_count()) where data is grouped by
   * vertex. For any vertex 'i' (where i ranges from 0 to m.vertex_count() - 1), all its adjacent
   * vertices are stored in a specific "block" within vtri:
   * - For vertex 'i', we access its neighbors in vtri from position:
   *    offset[i] to (offset[i] + degree[i] - 1).
   * - This block contains 'degree[i]' entries, representing the number of triangles that contain
   * vertex 'i'. This structure allows O(1) random access to the adjacency information of any
   * vertex.
   *
   * 2. DATA CONTENT (VERTEX NEIGHBORS)
   * Each entry in vtri[k] is a pair of indices {v_next, v_prev}.
   * - These two indices are the OTHER TWO vertices of the triangle that includes vertex 'i'.
   * - For example, if vertex 'i' is part of a triangle with vertices (i, A, B), the pair {A, B} is
   * stored.
   * - This is fundamental for FEM assembly, as these pairs identify the support of the basis
   * functions.
   */
}