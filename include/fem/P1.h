#pragma once

//#include "fem_matrix.h"
#include "sparse_matrix.h"
#include "mesh.h"

/* CSR variants */
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P);
void build_P1_mass_matrix(const Mesh &m, CSRMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, CSRMatrix &S);
void build_P1_mass_matrix(const Mesh &m, CSRPattern &P, CSRMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, CSRPattern &P, CSRMatrix &S);

/* SKL variants */
// void build_P1_SKLPattern(const Mesh &m, SKLPattern &P);
// void build_P1_mass_matrix(const Mesh &m, const SKLPattern &P, SKLMatrix &M);
// void build_P1_stiffness_matrix(const Mesh &m, const SKLPattern &P,
//			       SKLMatrix &S);

/* FEM matrix variants */
// void build_P1_mass_matrix(const Mesh &m, FEMatrix &M);
// void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S);
