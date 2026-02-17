#pragma once

#include "sparse_matrix.h"
#include "array.h"

/* Constructs the etree associated to a matrix A with sparsity pattern P.
Input: sparisty pattern P
vector parent of size n; should be same size as the matrix */
void construct_etree(CSRPattern &P, TArray<uint32_t> &parent);

void construct_L_sparsity_pattern(CSRPattern &A, CSRPattern &L, TArray<uint32_t> &parent);
