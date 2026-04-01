#if !defined(CHOLESKY_H)
#define CHOLESKY_H
#include "matrix.h"
#include "sparse_matrix.h"

using etree = TArray<uint32_t>;

//Phase 1: symbolic phase to determine the elimination tree
void symbolic_cholesky(const CSRMatrix &A, etree &T);

//Phase 2: symbolic factorization, compute the sparsity pattern of L, given the elimination tree
void L_pattern(const CSRMatrix &A,etree &T, CSRPattern &L_Pattern);

//3: numeric factorization, compute the numerical values of L given the sparsity pattern
void cholesky_fact(const CSRMatrix &A, CSRMatrix &L, const CSRPattern &L_Pattern);

//Phase 4: solve the linear system Ax = b using the Cholesky factorization, L
TArray<double> forward_sub(const CSRMatrix &L, const TArray<double> &b);
void backward_sub(const CSRMatrix &L, const TArray<double> &y, TArray<double> &x);
void cholesky_solve(const CSRMatrix &L, const TArray<double> &b, TArray<double> &x);

#endif /* CHOLESKY_H */