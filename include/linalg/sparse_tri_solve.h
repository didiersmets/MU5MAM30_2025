#include "sparse_matrix.h"

void sparse_tri_solve(CSRMatrix &L, double *__restrict b,double *__restrict x);
void sparse_tri_solve(CSRMatrix &L,
		double *b, uint32_t *b_ind, size_t b_size,
		double *x, uint32_t *x_ind, size_t x_size,
		uint32_t i);
