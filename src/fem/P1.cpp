#pragma once

#include <stdio.h>
#include <string.h>

// #include "fem_matrix.h"
// #include "sparse_matrix.h"

// File: src/fem/P1.cpp

// 1. Include own definitions (to check function prototypes)
#include "fem/P1.h"

// 2. Include dependencies for implementation:
#include "common/vec3.h"           // For Vec3, Vec3d, etc.
#include "fem/mass.h"              // For the local mass function mass()
#include "fem/stiffness.h"         // For the local stiffness function stiffness()
#include "matrix/fem_matrix.h"     // For the FEMatrix structure
#include "matrix/sparse_matrix.h"  // For CSRPattern, CSRMatrix, etc.
#include "mesh/mesh.h"             // For the Mesh structure

void build_P1_mass_matrix(const Mesh& m, FEMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, FEMatrix& S);

void build_P1_CSRPattern(const Mesh& m, CSRPattern& P);
void build_P1_mass_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& S);

void build_P1_SKLPattern(const Mesh& m, SKLPattern& P);
void build_P1_mass_matrix(const Mesh& m, const SKLPattern& P, SKLMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, const SKLPattern& P, SKLMatrix& S);

void build_P1_mass_matrix(const Mesh& m, FEMatrix& M)
{
  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();

  M.fem_type = FEMatrix::P1_cst;
  M.m        = &m;
  M.rows = M.cols = vtx_count;

  M.diag.resize(vtx_count);
  memset(M.diag.data, 0, vtx_count * sizeof(double));

  M.off_diag.resize(tri_count);
  const TArray<uint32_t>& idx = m.indices;
  for (size_t t = 0; t < tri_count; ++t)
  {
    uint32_t a  = idx[3 * t + 0];
    uint32_t b  = idx[3 * t + 1];
    uint32_t c  = idx[3 * t + 2];
    Vec3f    A  = m.positions[a];
    Vec3f    B  = m.positions[b];
    Vec3f    C  = m.positions[c];
    Vec3d    AB = {(double) B[0] - (double) A[0],
                   (double) B[1] - (double) A[1],
                   (double) B[2] - (double) A[2]};
    Vec3d    AC = {(double) C[0] - (double) A[0],
                   (double) C[1] - (double) A[1],
                   (double) C[2] - (double) A[2]};
    double   Mloc[2];
    mass(AB, AC, Mloc);
    // add the local contribution to the global matrix
    M.diag[a] += Mloc[0];
    M.diag[b] += Mloc[0];
    M.diag[c] += Mloc[0];

    M.off_diag[t] = Mloc[1];
  }
}