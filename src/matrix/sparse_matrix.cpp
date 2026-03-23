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
	uint32_t num_cols;		// num of nonzero cols in a row
	uint32_t index_in_data;	// index of the matrix entry in data array
	int j;					// column of the matrix entry and index input x
	for ( int i = 0; i < this->rows; i++ ) {
		y[i] = 0.0;
		num_cols = this->row_start[i+1] - row_start[i];
		for ( int j_index = 0; j_index < num_cols; j_index++ ) {
			index_in_data = this->row_start[i]+j_index;
			j = this->col[index_in_data];
			y[i] += this->data[index_in_data] * x[j];
		}
	}

	if ( symmetric ) {
		for ( int i = 0; i < this->rows; i++ ) {
			num_cols = this->row_start[i+1] - 1 - row_start[i];
			for ( int j_index = 0; j_index < num_cols; j_index++ ) {
				index_in_data = this->row_start[i]+j_index;
				j = this->col[index_in_data];
				y[j] += data[index_in_data] * x[i];
			}
		}
	}
	printf("mvp called\n");
}

double CSRMatrix::sum() const
{
	double res = 0.0;
	for (size_t k = 0; k < nnz; k++) {
		res += data[k];
	}
	// TODO: compress symmetric matrix such that only triangular matrix is stored
	if (symmetric) {
		res *= 2;
		for (size_t k = 0; k < rows; k++) {
			assert(col[row_start[k + 1] - 1] == k);
			res -= data[row_start[k + 1] - 1];
		}
	}
	return res;
}
