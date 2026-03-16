// include/linalg/sparse_cholesky.h
#pragma once

#include <stdint.h>
#include <stddef.h>
#include <vector>

#include "sparse_matrix.h"

struct ColAdjEntry {
    uint32_t row; // r
    size_t   k;   // index into L.data / L.col for entry (r, col)
};

void elimination_tree_from_csr_symmetric(const CSRPattern& P,
                                        std::vector<uint32_t>& parent);

void symbolic_cholesky_row_pattern(const CSRPattern& P,
                                   const std::vector<uint32_t>& parent,
                                   CSRPattern& PL);

void build_A_from_M_S(const CSRMatrix& M, const CSRMatrix& S,
                      double alpha, CSRMatrix& A);

void cholesky_factorize_rowwise(const CSRMatrix& A,
                                const CSRPattern& PL,
                                CSRMatrix& L);

void build_col_adjacency_lower(const CSRMatrix& L,
                               std::vector<std::vector<ColAdjEntry>>& col_adj);

void forward_solve_lower_csr(const CSRMatrix& L,
                             const double* b, double* y);

void backward_solve_lowerT_with_coladj(const CSRMatrix& L,
                                       const std::vector<std::vector<ColAdjEntry>>& col_adj,
                                       const double* y, double* x);