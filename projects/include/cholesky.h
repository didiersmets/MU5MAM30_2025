#include "sparse_matrix.h"
#include "array.h"

//assume L is prepared already with the correct sparsity patterna
void cholesky_factorization(CSRMatrix &A, CSRPattern&P,  CSRMatrix &L);

void solve_cholesky(CSRMatrix &L,  const double *__restrict b, double *__restrict x);
