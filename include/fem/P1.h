#pragma once

#include "fem_matrix.h"
#include "sparse_matrix.h"

//helper functions
bool find(const TArray<uint32_t> & arr, uint32_t target,  size_t start, size_t stop);


/* CSR variants */
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P);
void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P,
			       CSRMatrix &S);

/* FEM matrix variants */
void build_P1_mass_matrix(const Mesh &m, FEMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S);
