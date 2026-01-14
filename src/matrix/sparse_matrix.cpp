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
	// Initialisation du vecteur résultat
	for (size_t i = 0; i < rows; ++i) {
		y[i] = 0.0;
	}

	if (!symmetric) {
		// Cas général CSR
		for (size_t i = 0; i < rows; ++i) {
			size_t start = row_start[i];
			size_t stop  = row_start[i + 1];

			for (size_t k = start; k < stop; ++k) {
				y[i] += data[k] * x[col[k]];
			}
		}
	} else {
		// Cas symétrique
		for (size_t i = 0; i < rows; ++i) {
			size_t start = row_start[i];
			size_t stop  = row_start[i + 1];

			for (size_t k = start; k < stop; ++k) {
				uint32_t j = col[k];
				double a  = data[k];

				// Contribution diagonale
				if (i == j) {
					y[i] += a * x[i];
				}
				// Hors diagonale : contribution symétrique
				else {
					y[i] += a * x[j];
					y[j] += a * x[i];
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
