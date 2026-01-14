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
  /* (done) Your implementation goes here */

  for (uint32_t i=0; i<rows; ++i)
    y[i] = 0;

  for (uint32_t i=0; i<rows; ++i)
    for (uint32_t k=row_start[i]; k<row_start[i+1]; ++k) {
      uint32_t j = col[k];
      double a = data[k];

      y[i] += a * x[j];
      if ((symmetric) && (j != i))
	y[j] += a * x[i];
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
