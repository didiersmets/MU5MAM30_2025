#include "sparse_matrix.h"

#include <cassert>

/* CSRMatrix */

double& CSRMatrix::operator()(uint32_t i, uint32_t j)
{
  static double dummy = 0.0;
  assert(i < rows);
  size_t start = row_start[i];
  size_t stop  = row_start[i + 1];
  for (size_t k = start; k < stop; ++k)
  {
    if (col[k] == j)
      return data[k];
  }
  /* Not a valid matrix entry */
  assert(false);
  return dummy;
}

void CSRMatrix::mvp_full(const double* __restrict x, double* __restrict y) const
{
  // Safety check: this method expects a fully stored matrix
  assert(symmetric == false && "Error: mvp_full called on a symmetric-stored matrix.");

  // Initialize output vector to zero
  for (int i = 0; i < rows; i++)
    y[i] = 0.0;

  // Standard Sparse Matrix-Vector Multiplication (y = A * x)
  for (int i = 0; i < rows; i++)
  {
    double sum   = 0.0;
    size_t start = row_start[i];
    size_t stop  = row_start[i + 1];

    for (size_t k = start; k < stop; k++)
    {
      sum += data[k] * x[col[k]];
    }
    y[i] = sum;
  }
}

void CSRMatrix::mvp_symmetric(const double* __restrict x, double* __restrict y) const
{
  // Safety check: this method expects only the upper/lower triangular part
  assert(symmetric == true && "Error: mvp_symmetric called on a fully stored matrix.");

  // Initialize output vector to zero
  for (int i = 0; i < rows; i++)
    y[i] = 0.0;

  // Compute both A_upper and A_lower contributions in one pass , using y = A_upper * x and
  // y=A_upper^T * x
  for (int i = 0; i < rows; i++)
  {
    size_t start = row_start[i];
    size_t stop  = row_start[i + 1];

    for (size_t k = start; k < stop; k++)
    {
      uint32_t j     = col[k];
      double   value = data[k];

      // Standard contribution: y[i] += A[i][j] * x[j]
      y[i] += value * x[j];

      // Symmetric contribution: y[j] += A[j][i] * x[i]
      // We skip the diagonal (i == j) to avoid doubling the result
      if (i != j)
      {
        y[j] += value * x[i];
      }
    }
  }
}

double CSRMatrix::sum() const
{
  double res = 0.0;
  for (size_t k = 0; k < nnz; k++)
  {
    res += data[k];
  }
  if (symmetric)
  {
    res *= 2;
    for (size_t k = 0; k < rows; k++)
    {
      assert(col[row_start[k + 1] - 1] == k);
      res -= data[row_start[k + 1] - 1];
    }
  }
  return res;
}