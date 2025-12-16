#include "sparse_matrix.h"

/******************************************************************************
 *
 * Compressed Sparse Row matrix
 *
 *****************************************************************************/

// -----------------------------------
void CSRMatrix::mvp(const double *__restrict x, double *__restrict y) const
{
    /* Assume that y is not set to zero before */
    for (size_t i = 0; i < rows; i++)
    {
        y[i] = 0.;
        for (size_t k = row_start[i]; k < row_start[i + 1]; k++)
            y[i] += data[k] * x[col[k]];
    }
}

// -----------------------------------
double CSRMatrix::operator()(uint32_t i, uint32_t j)
{
    assert(i < rows);
    assert(j < cols);

    double val = 0.;

    if (row_start[i] < row_start[i + 1])
    {
        bool stop = false;
        size_t k = row_start[i];
        while (k < row_start[i + 1] && !stop)
        {
            if (j == col[k])
            {
                val = data[k];
                stop = true;
            }
            k++;
        }
    }
    return val;
}

// -----------------------------------
double CSRMatrix::sum() const
{
    double sum = 0;
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t k = row_start[i]; k < row_start[i + 1]; k++)
            sum += data[k];
    }

    return sum;
}
