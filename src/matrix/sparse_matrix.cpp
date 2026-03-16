#include "sparse_matrix.h"

/* CSRMatrix */

// double &CSRMatrix::operator()(uint32_t i, uint32_t j)
// {
// 	static double dummy = 0.0;
// 	assert(i < rows);
// 	size_t start = row_start[i];
// 	size_t stop = row_start[i + 1];
// 	for (size_t k = start; k < stop; ++k) {
// 		if (col[k] == j)
// 			return data[k];
// 	}
// 	/* Not a valid matrix entry */
// 	assert(false);
// 	return dummy;
// }

// void CSRMatrix::mvp(const double *__restrict x, double *__restrict y) const
// {
// 	/* Your implementation goes here */
// }

// double CSRMatrix::sum() const
// {
// 	double res = 0.0;
// 	for (size_t k = 0; k < nnz; k++) {
// 		res += data[k];
// 	}
// 	if (symmetric) {
// 		res *= 2;
// 		for (size_t k = 0; k < rows; k++) {
// 			assert(col[row_start[k + 1] - 1] == k);
// 			res -= data[row_start[k + 1] - 1];
// 		}
// 	}
// 	return res;
// }



#include "sparse_matrix.h"
#include <assert.h>

/* CSRMatrix 成员函数实现 */

// 1. 通过 (i, j) 坐标访问元素（主要用于调试，性能较低）
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
        double sum_val = 0.0;
        size_t start = row_start[i];
        size_t stop = row_start[i + 1];
        
        for (size_t k = start; k < stop; ++k) {
            sum_val += data[k] * x[col[k]];
        }
        y[i] = sum_val;
    }
}

double CSRMatrix::sum() const
{
    double res = 0.0;
    // 直接累加所有非零元素的值
    for (size_t k = 0; k < nnz; k++) {
        res += data[k];
    }

    // 如果矩阵是对称存储模式（只存了一半），才需要特殊处理
    // 但在 P1 模拟中，我们通常使用的是完整存储，所以直接返回 res
    if (symmetric) {
        res *= 2;
        /* Subtract diagonal entries once (they were counted twice above) */
        for (size_t k = 0; k < rows; k++) {
            /* Diagonal is stored as the last entry of each row */
            assert(col[row_start[k + 1] - 1] == k);
            res -= data[row_start[k + 1] - 1];
        }
    }
    return res;
}

/* 全局辅助函数 */

// 用于在 build 阶段查找 (i, j) 对应的 CSR 数据数组索引
size_t find_matrix_entry(const CSRMatrix &A, uint32_t i, uint32_t j) {
    uint32_t start = A.row_start[i];
    uint32_t end   = A.row_start[i + 1];
    
    for (uint32_t k = start; k < end; ++k) {
        if (A.col[k] == j) {
            return (size_t)k;
        }
    }
    // 如果执行到这里，说明你的 CSRPattern 漏掉了某个非零项
    assert(false && "Entry not found in CSR pattern!");
    return 0;
}