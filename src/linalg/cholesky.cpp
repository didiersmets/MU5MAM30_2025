#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>
#include "tiny_blas.h"
#include "cholesky.h"
#include <vector>

void symbolic_cholesky(const CSRMatrix &A, etree &T)
{
    size_t n = A.rows;
    constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    etree ancestors(n, INVALID);
    T.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        T[i] = INVALID;
        uint32_t start = A.row_start[i];
        uint32_t stop = A.row_start[i + 1];

        for (size_t k = start; k < stop; ++k) {
            uint32_t j = A.col[k];
            if (j < i) {
                uint32_t jroot = j;
                while (jroot != INVALID && ancestors[jroot] != (uint32_t)i) {
                    uint32_t next = ancestors[jroot];
                    ancestors[jroot] = (uint32_t)i;
                    if (next == INVALID) {
                        T[jroot] = (uint32_t)i;
                    }
                    jroot = next;
                }
            }
        }
    }
}

void L_pattern(const CSRMatrix &A, etree &T, CSRPattern &L_Pattern)
{
    size_t n = A.rows;
    constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    L_Pattern.rows = n;
    L_Pattern.cols = n;
    L_Pattern.row_start.resize(n + 1);
    TArray<TArray<uint32_t>> L_cols(n);
    TArray<uint32_t> mark(n);
    for (size_t i = 0; i < n; ++i){
        mark[i] = i;
        uint32_t start = A.row_start[i];
        uint32_t stop = A.row_start[i+1];
        for (size_t k = start; k < stop; ++k){
            if (A.col[k] < i){
                uint32_t j = A.col[k];
                while(j != INVALID && mark[j] != i){
                    L_cols[i].push_back(j);
                    mark[j] = i;
                    j = T[j];
                }
            }
        }
        L_cols[i].push_back(i);
    }

    uint32_t nnz = 0;
    for (size_t i = 0; i < n; ++i){
        L_Pattern.row_start[i] = nnz;
        nnz += L_cols[i].size;
    }
    L_Pattern.nnz = nnz;
    L_Pattern.row_start[n] = nnz;
    L_Pattern.col.resize(nnz);
    for (size_t i = 0; i < n; ++i){
        uint32_t start = L_Pattern.row_start[i];
        for (size_t k = 0; k < L_cols[i].size; ++k){
            L_Pattern.col[start + k] = L_cols[i][k];
        }
    }
}

void cholesky_fact(const CSRMatrix &A, CSRMatrix &L, const CSRPattern &L_Pattern) {
    size_t n = A.rows;
    L.rows = n;
    L.cols = n;
    L.nnz = L_Pattern.nnz;

    L.row_start = L_Pattern.row_start.data;
    L.col = L_Pattern.col.data;
    L.data.resize(L.nnz);

    std::vector<double> work(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t k = A.row_start[i]; k < A.row_start[i + 1]; ++k) {
            if (A.col[k] <= i) {
                work[A.col[k]] = A.data[k];
            }
        }

        size_t l_start = L.row_start[i];
        size_t l_stop = L.row_start[i + 1];

        for (size_t k = l_start; k < l_stop; ++k) {
            uint32_t j = L.col[k];

            if (j < i) {
                double l_jj = 1.0;
                size_t m_start = L.row_start[j];
                size_t m_stop = L.row_start[j + 1];

                for (size_t m = m_start; m < m_stop; ++m) {
                    uint32_t col_m = L.col[m];
                    if (col_m < j) {
                        work[j] -= work[col_m] * L.data[m];
                    } else if (col_m == j) {
                        l_jj = L.data[m];
                    }
                }
                
                work[j] /= l_jj;
                L.data[k] = work[j];

            } else if (j == i) {
                double sum_sq = 0.0;
                for (size_t m = l_start; m < k; ++m) {
                    sum_sq += L.data[m] * L.data[m];
                }
                L.data[k] = std::sqrt(work[i] - sum_sq);
            }
        }

        for (size_t k = l_start; k < l_stop; ++k) {
            work[L.col[k]] = 0.0;
        }
    }
}

//Maintenant on va implémenter la forward et la backward substitution

TArray<double> forward_sub(const CSRMatrix &L, const TArray<double> &b){
    size_t n = L.rows;
    TArray<double> y(n, 0.0);

    for(size_t i = 0; i < n; ++i){
        double sum = 0.0;
        size_t start = L.row_start[i];
        size_t stop = L.row_start[i + 1];
        for(size_t k = start; k < stop; ++k){
            uint32_t j = L.col[k];
            if(j < i){
                sum += L.data[k] * y[j];
            }
            else if(j == i){
                y[i] = b[i] - sum;
                y[i] /= L.data[k];
                break;
            }
            }
        }
    return y;
}


void backward_sub(const CSRMatrix &L, const TArray<double> &y, TArray<double> &x){
    size_t n = L.rows;
    x = y; // copy y to x

    for (size_t i = n; i-- > 0; ) {
        //Division par la diagonale qui est le dernier élément de la ligne i
        size_t diag_k = L.row_start[i + 1] - 1;
        x[i] /= L.data[diag_k];

        //On distribue la contibution de x[i] sur tous les x[j] (avec j < i)
        size_t start = L.row_start[i];
        for (size_t k = start; k < diag_k; ++k) {
            uint32_t j = L.col[k];
            x[j] -= L.data[k] * x[i];
        }
    }
}

void cholesky_solve(const CSRMatrix &L, const TArray<double> &b, TArray<double> &x){
    backward_sub(L, forward_sub(L,b), x);
}

