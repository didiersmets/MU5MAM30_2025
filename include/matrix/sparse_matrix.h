#pragma once

#include "common/array.h"
#include "matrix.h"

#include <stdint.h>

/******************************************************************************
 *
 * Compressed Sparse Row matrix
 *
 *****************************************************************************/

struct CSRPattern
{
  bool   symmetric;
  size_t rows;
  size_t cols;
  size_t nnz;
  /* Non zero entries on line i (0 <= i < rows)
   * are stored at indices row_start(i) <= k < row_start(i + 1).
   * Corresponding column indices are read into col(k). */
  // offset for every row: row_start[i+1] - row_start[i] = number of non zero entries in row i
  TArray<uint32_t> row_start; /* Size = nrows + 1 */
  // indices of columns for every non zero entry
  TArray<uint32_t> col; /* Size = nnz */
};

struct CSRMatrix : public Matrix
{
  bool           symmetric = false;
  size_t         nnz;       /* Number of (non zero) entries */
  uint32_t*      row_start; /* pointer to the corresponding data in pattern */
  uint32_t*      col;       /* pointer to the corresponding data in pattern */
  TArray<double> data;      /* Size = nnz  */
  void           mvp(const double* __restrict x, double* __restrict y) const;
  double         sum() const;
  double&        operator()(uint32_t i, uint32_t j);
};
