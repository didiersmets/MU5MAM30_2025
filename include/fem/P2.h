//
// Created by aurel on 11/02/2026.
//
#include "mesh.h"
#include "sparse_matrix.h"
#include "fem_matrix.h"

#ifndef MU5MAM30_P2_H
#define MU5MAM30_P2_H

#include "adjacency_P2.h"

/* CSR variants for P2 */
// NOUVEAU : On ajoute EdgeAdjacency
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P, const EdgeAdjacency &edges);
void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M, const EdgeAdjacency &edges);
void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S, const EdgeAdjacency &edges);

/* FEM matrix variants for P2 */
// (On les garde déclarées au cas où tu passerais USE_FEM_MATRIX à true un jour)
void build_P2_mass_matrix(const Mesh &m, FEMatrix &M);
void build_P2_stiffness_matrix(const Mesh &m, FEMatrix &S);

#endif //MU5MAM30_P2_H