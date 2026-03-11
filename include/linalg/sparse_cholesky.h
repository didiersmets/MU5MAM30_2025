#include "array.h"
#include "sparse_matrix.h"

void elimination_tree(const CSRMatrix &A,TArray<uint32_t> &parent);
void cholesky_sparsity_pattern(const CSRMatrix &A,const TArray<uint32_t> &parent,CSRPattern &cho_patt);
