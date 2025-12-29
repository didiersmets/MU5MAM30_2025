#include "sparse_matrix.h"

/* CSRMatrix */

double &CSRMatrix::operator()(uint32_t i, uint32_t j)
{
	static double dummy = 0.0;
	assert(i < rows);
	size_t start = row_start[i];
	size_t stop = row_start[i + 1];
	for (size_t k = start; k < stop; ++k) {
		if (col[k] == j)
			return data[k];
	}
	/* Not a valid matrix entry */
	assert(false);
	return dummy;
}

void CSRMatrix::mvp(const double *__restrict x, double *__restrict y) const
{
	/* Your implementation goes here */
	/* 
	Performs sparse matrix-vector multiplication (y = A * x) using CSR format. 
	For each row, it accumulates contributions from non-zero entries: 
	y[i] += data[k] * x[col[k]]. 
	If the matrix is symmetric (upper triangular stored), it adds mirrored contributions 
	for off-diagonal elements: y[col[k]] += data[k] * x[i]. 
 	*/

	for( int i = 0; i < rows; i++ ){
		y[i] = 0.0;
		size_t start = row_start[i];
		size_t stop  = row_start[i + 1];
		for( size_t k = start; k < stop; k++ ){
			y[i] += data[k] * x[ col[k] ];
		}
	}
	if( symmetric ){
		//we need to do the equivalent of A^T * x
		for( int i = 0; i < rows; i++ ){
			size_t start = row_start[i];
			size_t stop  = row_start[i + 1];
			for( size_t k = start; k < stop - 1; k++ ){ //-1 to avoid the diagonal
				y[ col[k] ] += data[k] * x[i];
			}
		}
	}
}

double CSRMatrix::sum() const
{
	double res = 0.0;
	for (size_t k = 0; k < nnz; k++) {
		res += data[k];
	}
	if (symmetric) {
		res *= 2;
		for (size_t k = 0; k < rows; k++) {
			assert(col[row_start[k + 1] - 1] == k);
			res -= data[row_start[k + 1] - 1];
		}
	}
	return res;
}
