#pragma once

#include "fem/P1.h"

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

/*
void build_P1_SKLPattern(const Mesh& m, SKLPattern& P);
void build_P1_mass_matrix(const Mesh& m, const SKLPattern& P, SKLMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, const SKLPattern& P, SKLMatrix& S);
*/

// global Mass and Stiffness matrix builders using CSR patterns

void build_P1_CSRPattern(const Mesh& m, CSRPattern& P)
{
  /* Your implementation goes here.
   * Use a VTAdjacency structure (see src/mesh/adjacency.cpp) to help build the pattern.
   */

  // CSR Pattern of a matrix to develop the pattern (cols indices, row offset)
}

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

  /* Your implementation goes here */
  // assemble local matrix M_loc and then add in the global matrix M
}

void build_P1_stiffness_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& S)
{
  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();
  assert(P.row_start.size == vtx_count + 1);

  S.symmetric = true;
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
