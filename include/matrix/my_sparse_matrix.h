#pragma once

#include <stdint.h>
#include "array.h"
struct MyCSRPattern {
	TArray<uint32_t> row_start;
	TArray<uint32_t> col;
};
struct MyCSRMatrix {

	size_t rows;
	size_t cols;
	bool symmetric;
	size_t nnz; 

	uint32_t *row_start;
	uint32_t *col;
	TArray<double> data; 
	void mvp(const double *__restrict x, double *__restrict y) const;
	double sum() const;
	double &operator()(uint32_t i, uint32_t j);
};

/* Compatibility aliases so code using the project's CSR types or the
 * custom MyCSR types both compile. */
typedef MyCSRPattern CSRPattern;
typedef MyCSRMatrix CSRMatrix;