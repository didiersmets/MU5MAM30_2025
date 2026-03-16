// matrix_sum.h
#pragma once

#include "matrix.h"
#include "tiny_blas.h"
#include <vector>

struct MatrixSum : public Matrix {
    const Matrix &A;
    const Matrix &B;
    double alpha;  // coefficient for B

    MatrixSum(const Matrix &A, const Matrix &B, double alpha, size_t N)
        : A(A), B(B), alpha(alpha)
    {
        rows = N;
        cols = N;
    }

    void mvp(const double *__restrict x, double *__restrict y) const override
    {
        std::vector<double> Bx(rows, 0.0);
        A.mvp(x, y);                          // y = A*x
        B.mvp(x, Bx.data());                  // Bx = B*x
        blas_axpy(alpha, Bx.data(), y, rows);  // y += alpha*B*x
    }

    double sum() const override { return 0.0; }
};