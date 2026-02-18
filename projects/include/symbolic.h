#pragma once

#include "sparse_matrix.h"
#include "array.h"

/* Constructs the etree associated to a matrix A with sparsity pattern P.
Input: P - Sparsity pattern of SPD matrix A
parent - vector of size n; should be same size as the matrix where the etree will be stored*/
void construct_etree(CSRPattern &P, TArray<uint32_t> &parent);

/* Construct the sparsity pattern of Cholesky factor L of SPD matrix A
Input: A - sparsity pattern of A, matrix to be factorized
L - sparsity pattern of L (to be constructed)
parent - vector describing the etree associated to A
*/
void construct_L_sparsity_pattern(CSRPattern &A, CSRPattern &L, TArray<uint32_t> &parent);
