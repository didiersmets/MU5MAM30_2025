#include "sparse_matrix.h"
#include "math_utils.h"

void CSRMatrix::mvp(const double *__restrict x, double *__restrict y) const
{
    for(size_t i = 0; i < rows; i++){
        double sum=0.0;
        for(size_t tmp= row_start[i]; tmp < row_start[i+1]; tmp++){
            uint32_t j = col[tmp];
            sum += data[tmp]* x[j];
        }

        y[i]= sum;
    }
    
}

double CSRMatrix::sum() const
{
    double sum = 0.0;
    for(size_t i=0; i< nnz; i++){
        sum += data[i];
    }

    return sum;
}

double &CSRMatrix::operator()(uint32_t i, uint32_t j)
{
    ASSERT(i < rows);
    ASSERT(j < cols);

    for (uint32_t k = row_start[i]; k < row_start[i+1]; k++) {
        if (col[k] == j) {
            return data[k];
        }
    }

    ASSERT_ALWAYS(false && "CSRMatrix: accessing entry not in sparsity pattern.");
     
    static double dummy = 0.0;
    return dummy;
}
