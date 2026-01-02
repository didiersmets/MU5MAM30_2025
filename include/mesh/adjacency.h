#include "common/array.h"
#include "mesh/mesh.h"

#include <stdint.h>

/*****************************************************************************
 * The mesh structure (mesh.h) that we use does not provide an immediate way
 * to determine which triangles are incident to a given vertex.
 * Since that information is helpful in some situations (we shall use it e.g.
 * to determine the CSR sparsity pattern of the associated mass and stiffness
 * matrices), we provide a dedicated structure, called VTAdjacency, than can
 * be initialized on the fly from a given existing mesh, and which encodes that
 * information.
 *
 * The incident triangles to vertex number "a" are not encoded by their
 * in the triangles array, but instead by listing the indices of the two
 * remaining vertices, in order : more precisely, if one such triangle is
 * "abc" then the VTri struct below will record the pair (b,c), because "a"
 * is known (no need to record it again), "b" is the "next" vertex and "c" is
 * the "prev" vertex in the triangle orientation.
 *
 * Once loaded, the 3 arrays containing information are :
 *
 * - degree : an array of size m.vertex_count(), degree[a] is the number of
 *            triangles incident to vertex a.
 * - offset : an array of size m.vertex_count(), this is just an accumulation
 *            of the degrees or earlier vertices, it is used to inform in which
 *            starting position (in the array vtri below) the adjacency info
 *            of a given vertex can be found. In simple math terms :
 *            offset[0] = 0 and offset[k] = \sum_{i = 0}^{k-1} degree[i]
 * - vtri   : an array of size 3 * m.triangle_count() [indeed each triangle
 *            lists 3 vertices so the total sum of the degrees is 3 times as
 *            much as the number of vertices]. The VTri structs (as explained
 *            above) related to vertex k are found in the vtri array at
 *            positions j such that offset[k] <= j < offser[k] + degree[k].
 *
 *****************************************************************************/

/* Vertex to Triangle adjacency table */
struct VTAdjacency
{
  struct VTri
  {
    uint32_t next;
    uint32_t prev;
  };
  TArray<uint32_t> degree;
  TArray<uint32_t> offset;
  TArray<VTri>     vtri;
  /* The constructor that you need to implement in src/adjacency.cpp */
  VTAdjacency(const Mesh& m);
};