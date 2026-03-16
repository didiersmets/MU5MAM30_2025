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
    for (size_t i = 0; i < rows; ++i) {
        y[i] = 0.0;
    }

    for (size_t i = 0; i < rows; ++i) {
        const size_t begin = row_start[i];
        const size_t end   = row_start[i + 1];

        for (size_t k = begin; k < end; ++k) {
            assert(k < nnz);
            assert(col[k] < cols);

            y[i] += data[k] * x[col[k]];
        }
    }

    if (symmetric) {
        for (size_t i = 0; i < rows; ++i) {
            const size_t begin = row_start[i];
            const size_t end   = row_start[i + 1];

            for (size_t k = begin; k < end; ++k) {
                const size_t j = col[k];

                if (j != i) {
                    y[j] += data[k] * x[i];
                }
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
