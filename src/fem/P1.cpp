#include "fem/P1.h"

#include "fem/mass.h"
#include "fem/stiffness.h"
#include "math.h"
#include "mesh/adjacency.h"
#include "mesh/mesh.h"

#include <algorithm>
#include <stdio.h>
#include <string.h>

/*
#include "P1.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness.h"
*/

void build_P1_CSRPattern(const Mesh& m, CSRPattern& P);
void build_P1_mass_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& S);

/**
 * Helper function: linear search to avoid duplicate column indices in the same row.
 * Returns true if the index 'x' is already present in the current row segment.
 */
static bool find(uint32_t x, uint32_t* start, size_t count)
{
  for (size_t i = 0; i < count; ++i)
  {
    if (start[i] == x)
      return true;
  }
  return false;
}

/**
 * Constructs the Compressed Sparse Row (CSR) pattern for P1 finite elements.
 * Each row 'a' represents vertex 'a' and its connections to other vertices.
 */
void build_P1_CSRPattern(const Mesh& m, CSRPattern& P)
{
  /*
    build the sparse matrix pattern to store information on how
    all the vertixes are connected.
    The main idea is to use the Adjency pattern built before:

    for the row array we can just use the offset array we computed
    in the constructor of Adjency since it already tells us on
    each row how many entries we have, since in each row we will have
    as many entries as the degree of that specific vtx.

    for the column array we instead need to loop over the vtri array
    in this array we will have informations on which are the vtxs to which
    our vtx is connected to, for each vtx we need to loop from
    vtri[ offset[vtx] ] to vtri[ offset[vtx] + degree ]

   */
  P.symmetric      = true;
  P.rows           = m.vertex_count();
  size_t vtx_count = m.vertex_count();
  P.cols           = m.vertex_count();
  VTAdjacency adj(m);

  P.row_start.resize(P.rows + 1);
  // we start with the worst case scenario where we fill everything
  P.col.resize(P.cols + 3 * m.triangle_count());

  int nnz = 0;
  for (int a = 0; a < P.rows; a++)
  {
    // also try with offset array in adjacency
    //  P.row_start[a] = adj.offset[a]
    P.row_start[a] = nnz;

    int init_nnz = nnz;

    int start = adj.offset[a];
    int stop  = start + adj.degree[a];

    for (int i = start; i < stop; i++)
    {
      // retrieve the entries from the vtri
      int b = adj.vtri[i].next;
      int c = adj.vtri[i].prev;

      /*
      since we are iterating for the number of a connected to vtx
      and not on the triangles that vtx is part of we will end up with
      redundant entries in the column array, to avoid that we need to see
      if the entries are already present
      */
      bool b_present    = false;
      bool c_present    = false;
      int  tmp_add_elem = nnz - init_nnz;

      for (int j = 0; j < tmp_add_elem; j++)
      {
        if (P.col[init_nnz + j] == b)
          b_present = true;
        if (P.col[init_nnz + j] == c)
          c_present = true;
      }
      // since the mtx is symmetric we only add the lower diagonal elements
      // so we check b < a and c < a
      if (b < a && !b_present)
        P.col[nnz++] = b;
      if (c < a && !c_present)
        P.col[nnz++] = c;
    }
    P.col[nnz++] = a;  // add the diagonal entry
  }

  P.col.resize(nnz);

  // Finalize row_start with the actual total number of non-zero elements
  P.row_start[vtx_count] = nnz;
  P.col.resize(nnz);
  P.col.shrink_to_fit();

  /* PHASE 2: SORT COLUMN INDICES */
  /* Standard CSR format requires column indices to be sorted for each row. */
  for (size_t a = 0; a < vtx_count; ++a)
  {
    uint32_t* __restrict to_sort = &P.col[P.row_start[a]];
    size_t count                 = P.row_start[a + 1] - P.row_start[a];

    /* Insertion sort: very efficient for small arrays (vertex degrees are usually low) */
    for (size_t k = 1; k < count; ++k)
    {
      size_t j = k;
      while (j && (to_sort[j - 1] > to_sort[j]))
      {
        uint32_t tmp   = to_sort[j - 1];
        to_sort[j - 1] = to_sort[j];
        to_sort[j]     = tmp;
        j--;
      }
    }
  }
}

/*
struct CSRPattern
{
  bool   symmetric;
  size_t rows;
  size_t cols;
  size_t nnz;
   Non zero entries on line i (0 <= i < rows)
   * are stored at indices row_start(i) <= k < row_start(i + 1).
   * Corresponding column indices are read into col(k).
  // offset for every row: row_start[i+1] - row_start[i] = number of non zero entries in row i
  TArray<uint32_t> row_start;    Size = nrows + 1
  // indices of columns for every non zero entry
  TArray<uint32_t> col;       Size = nnz
  }
  */

void build_P1_mass_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& M)
{
  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();
  assert(P.row_start.size == vtx_count + 1);

  M.symmetric = true;
  M.rows = M.cols = vtx_count;
  M.nnz           = P.col.size;
  M.row_start     = P.row_start.data;
  M.col           = P.col.data;
  M.data.resize(M.nnz);
  for (size_t i = 0; i < M.nnz; ++i)
  {
    M.data[i] = 0.0;
  }

  // global indices of the mesh triangles
  const TArray<uint32_t>& idx = m.indices;

  // Main assembly loop: iterate over each triangle (element) in the mesh
  for (size_t t = 0; t < tri_count; ++t)
  {
    // Retrieve global indices for the three vertices of triangle 't'
    uint32_t v[3] = {idx[3 * t + 0], idx[3 * t + 1], idx[3 * t + 2]};

    // Retrieve the 3D positions of the triangle's vertices
    Vec3f A_pos = m.positions[v[0]];
    Vec3f B_pos = m.positions[v[1]];
    Vec3f C_pos = m.positions[v[2]];

    // Edge vectors
    Vec3d AB = {(double) B_pos[0] - A_pos[0],
                (double) B_pos[1] - A_pos[1],
                (double) B_pos[2] - A_pos[2]};
    Vec3d AC = {(double) C_pos[0] - A_pos[0],
                (double) C_pos[1] - A_pos[1],
                (double) C_pos[2] - A_pos[2]};

    // compute the local mass matrix for the triangle ABC
    double Mloc[2];
    mass(AB, AC, Mloc);

    // GLOBAL ASSEMBLY : add the local mass matrix Mloc into the global matrix M
    // Iterate over all pairs (i, j) of the 3x3 local element matrix
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        // contribution to global matrix M at position (r, c)

        uint32_t r = v[i];  // Global row index
        uint32_t c = v[j];  // Global column index

        /* * Since the matrix is symmetric and we use Lower Triangular storage,
         * we only process entries where row index 'r' >= column index 'c'.
         * Entries where r < c belong to the upper triangle and are omitted.
         */
        if (r >= c)
        {
          // Determine value: use diagonal term if r==c, else use off-diagonal term
          double val_to_add = (r == c) ? Mloc[0] : Mloc[1];

          // Find the position 'k' in the M.data array corresponding to column 'c'
          // We only search within the range defined for row 'r'
          uint32_t start = M.row_start[r];
          uint32_t end   = M.row_start[r + 1];

          for (uint32_t k = start; k < end; ++k)
          {
            if (M.col[k] == c)
            {
              // accumulate the local contribution
              M.data[k] += val_to_add;
              break;
            }
          }
        }
      }
    }
  }
}

void build_P1_stiffness_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& S)
{
  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();

  // Ensure the pattern matches the mesh size
  assert(P.row_start.size == vtx_count + 1);

  // 1. Initialize CSR Matrix properties
  // We only store the lower triangular part (symmetric = true)
  S.symmetric = true;
  S.rows      = vtx_count;
  S.cols      = vtx_count;
  S.nnz       = P.col.size;

  // 2. Bind the matrix to the topology provided by the CSRPattern
  S.row_start = P.row_start.data;
  S.col       = P.col.data;

  // 3. Allocate and zero-initialize the numerical data array
  S.data.resize(S.nnz);
  for (size_t i = 0; i < S.nnz; ++i)
  {
    S.data[i] = 0.0;
  }

  /* ASSEMBLY PHASE */
  double Sloc[6];  // Local stiffness matrix components: [Saa, Sbb, Scc, Sab, Sbc, Sca]
  const TArray<uint32_t>& idx = m.indices;

  for (size_t t = 0; t < tri_count; ++t)
  {
    // Retrieve global indices for the three vertices of triangle 't'
    uint32_t v[3] = {idx[3 * t + 0], idx[3 * t + 1], idx[3 * t + 2]};

    // Retrieve the 3D positions of the triangle's vertices
    Vec3f A_pos = m.positions[v[0]];
    Vec3f B_pos = m.positions[v[1]];
    Vec3f C_pos = m.positions[v[2]];

    // Edge vectors
    Vec3d AB = {(double) B_pos[0] - A_pos[0],
                (double) B_pos[1] - A_pos[1],
                (double) B_pos[2] - A_pos[2]};
    Vec3d AC = {(double) C_pos[0] - A_pos[0],
                (double) C_pos[1] - A_pos[1],
                (double) C_pos[2] - A_pos[2]};

    // Calculate local stiffness contributions for triangle ABC
    stiffness(AB, AC, Sloc);

    /* * MAPPING TO LOWER TRIANGULAR STORAGE:
     * To respect the 'S.symmetric = true' constraint, we must ensure the row index (r)
     * is always >= the column index (c). We achieve this by taking the max and min
     * of each global vertex pair.
     */
    struct Interaction
    {
      uint32_t row, col;
      double   value;
    };
    Interaction pairs[6] = {
      {v[0], v[0], Sloc[0]},                                  // Saa (Diagonal)
      {v[1], v[1], Sloc[1]},                                  // Sbb (Diagonal)
      {v[2], v[2], Sloc[2]},                                  // Scc (Diagonal)
      {std::max(v[0], v[1]), std::min(v[0], v[1]), Sloc[3]},  // Interaction A-B
      {std::max(v[1], v[2]), std::min(v[1], v[2]), Sloc[4]},  // Interaction B-C
      {std::max(v[2], v[0]), std::min(v[2], v[0]), Sloc[5]}   // Interaction C-A
    };

    // Scatter each of the 6 local interactions into the global CSR data array
    for (int k_loc = 0; k_loc < 6; ++k_loc)
    {
      uint32_t r   = pairs[k_loc].row;
      uint32_t c   = pairs[k_loc].col;
      double   val = pairs[k_loc].value;

      // Locate the entry for column 'c' within the compressed row 'r'
      uint32_t row_begin = S.row_start[r];
      uint32_t row_end   = S.row_start[r + 1];

      for (uint32_t k = row_begin; k < row_end; ++k)
      {
        if (S.col[k] == c)
        {
          // Atomically add the local triangle contribution to the global matrix entry
          S.data[k] += val;
          break;  // Found the column, move to the next interaction
        }
      }
    }
  }
}
