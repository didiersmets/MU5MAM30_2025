#pragma once

#include "fem/P1.h"

#include "fem/mass.h"
#include "mesh/adjacency.h"
#include "mesh/mesh.h"

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

void build_P1_mass_matrix(const Mesh& m, FEMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, FEMatrix& S);

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
  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();

  // The matrix represents interactions between mesh vertices (V x V).
  // We store only the lower triangle (including diagonal) to save memory.
  P.symmetric = true;
  P.rows = P.cols = vtx_count;
  P.row_start.resize(vtx_count + 1);

  // Get vertex-to-triangle adjacency data (provides triangles sharing each vertex)
  VTAdjacency adj(m);

  /* Initial memory allocation estimate:
   * Each triangle (tri_count) has 3 edges, plus the diagonal for each vertex.
   */
  size_t max_nnz = 3 * tri_count + vtx_count;
  P.col.resize(max_nnz);

  /* PHASE 1: FILL THE PATTERN (UNORDERED) */
  size_t nnz = 0;
  for (size_t a = 0; a < vtx_count; ++a)
  {
    // P.row_start[a] stores the offset (starting position) of row 'a' in the flat P.col array.
    P.row_start[a] = nnz;

    // Pointer to the beginning of row 'a' in P.col for quick searching
    uint32_t* start   = &P.col[nnz];
    size_t    nnz_loc = 0;  // Local counter for elements added to this specific row

    // Access the block of triangles associated with vertex 'a' from VTAdjacency
    uint32_t kstart = adj.offset[a];
    uint32_t kstop  = kstart + adj.degree[a];

    for (size_t k = kstart; k < kstop; ++k)
    {
      // b and c are the neighbors of 'a' in the current triangle k.
      uint32_t b = adj.vtri[k].next;
      uint32_t c = adj.vtri[k].prev;

      /* * Check neighbor 'b':
       * 1. 'b < a' ensures we only store the lower triangular part of the matrix.
       * 2. '!find' ensures we don't add the same neighbor twice (shared edges).
       */
      if (b < a && !find(b, start, nnz_loc))
      {
        P.col[nnz++] = b;
        nnz_loc++;
      }
      /* Check neighbor 'c' with the same logic */
      if (c < a && !find(c, start, nnz_loc))
      {
        P.col[nnz++] = c;
        nnz_loc++;
      }
    }
    // Add vertex 'a' itself (diagonal element).
    P.col[nnz++] = a;
  }

  // Finalize row_start with the total number of non-zero elements
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
  /* Non zero entries on line i (0 <= i < rows)
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
  assert(P.row_start.size == vtx_count + 1);

  S.symmetric = true;  // saving only the lower triangular part
  S.rows = S.cols = vtx_count;
  S.nnz           = P.col.size;
  S.row_start     = P.row_start.data;
  S.col           = P.col.data;
  S.data.resize(S.nnz);
  for (size_t i = 0; i < S.nnz; ++i)
  {
    S.data[i] = 0.0;
  }

  /* Your implementation goes here */
  // assemble local matrix S_loc and then add in the global matrix S
}
