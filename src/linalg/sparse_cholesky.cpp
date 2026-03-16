// src/linalg/sparse_cholesky.cpp
#include "sparse_cholesky.h"

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <limits>

static constexpr uint32_t NIL = std::numeric_limits<uint32_t>::max();

static inline size_t find_pos_in_row_binary(const CSRMatrix& A, uint32_t i, uint32_t j)
{
    size_t lo = (size_t)A.row_start[i];
    size_t hi = (size_t)A.row_start[i + 1];
    while (lo < hi) {
        size_t mid = (lo + hi) >> 1;
        uint32_t c = A.col[mid];
        if (c < j) lo = mid + 1;
        else       hi = mid;
    }
    assert(lo < (size_t)A.row_start[i + 1] && A.col[lo] == j);
    return lo;
}

// dot over k < j of L(i,k) * L(j,k)
static inline double dot_row_row_until_j(const CSRMatrix& L, uint32_t i, uint32_t j)
{
    size_t pi   = (size_t)L.row_start[i];
    size_t pj   = (size_t)L.row_start[j];
    size_t endi = (size_t)L.row_start[i + 1] - 1; // exclude diag (last)
    size_t endj = (size_t)L.row_start[j + 1] - 1; // exclude diag

    double s = 0.0;

    while (pi < endi && pj < endj) {
        uint32_t ci = L.col[pi];
        uint32_t cj = L.col[pj];

        // only k < j
        if (ci >= j || cj >= j) break;

        if (ci == cj) {
            s += L.data[pi] * L.data[pj];
            ++pi; ++pj;
        } else if (ci < cj) {
            ++pi;
        } else {
            ++pj;
        }
    }
    return s;
}

/* Algorithm 4.2 */

void elimination_tree_from_csr_symmetric(const CSRPattern& P,
                                        std::vector<uint32_t>& parent)
{
    const uint32_t n = (uint32_t)P.rows;
    parent.assign(n, NIL);

    std::vector<uint32_t> ancestor(n, NIL);

    for (uint32_t i = 0; i < n; ++i) {
        // below diagonal entries in row i are all cols except the last (diag=i)
        size_t start = (size_t)P.row_start[i];
        size_t stop  = (size_t)P.row_start[i + 1];
        if (stop == start) continue;

        // assume diag is present and sorted so diag is last
        stop -= 1;

        for (size_t kk = start; kk < stop; ++kk) {
            uint32_t j = P.col[kk];
            assert(j < i);

            uint32_t jroot = j;
            while (ancestor[jroot] != NIL && ancestor[jroot] != i) {
                uint32_t l = ancestor[jroot];
                ancestor[jroot] = i;   // path compression
                jroot = l;
            }
            if (ancestor[jroot] == NIL) {
                ancestor[jroot] = i;
                parent[jroot] = i;
            }
        }
    }
}

/* Algorithm 4.3 -> build CSRPattern of L */

void symbolic_cholesky_row_pattern(const CSRPattern& P,
                                   const std::vector<uint32_t>& parent,
                                   CSRPattern& PL)
{
    const uint32_t n = (uint32_t)P.rows;
    assert(parent.size() == n);

    std::vector<std::vector<uint32_t>> rowL(n);
    std::vector<int32_t> mark(n, -1);

    for (uint32_t i = 0; i < n; ++i) {
        rowL[i].clear();
        mark[i] = (int32_t)i;

        size_t start = (size_t)P.row_start[i];
        size_t stop  = (size_t)P.row_start[i + 1];
        if (stop == start) {
            continue;
        }
        stop -= 1;
        for (size_t kk = start; kk < stop; ++kk) {
            uint32_t k = P.col[kk];
            assert(k < i);

            uint32_t j = k;
            while (j != NIL && mark[j] != (int32_t)i) {
                mark[j] = (int32_t)i;
                rowL[i].push_back(j);
                j = parent[j];
            }
        }

        std::sort(rowL[i].begin(), rowL[i].end());
    }

    // Build PL (row-wise, with diag at the end of each row)
    PL.symmetric = true;
    PL.rows = PL.cols = n;

    PL.row_start.resize((size_t)n + 1);
    PL.row_start[0] = 0;

    size_t nnz = 0;
    for (uint32_t i = 0; i < n; ++i) {
        nnz += rowL[i].size() + 1; // + diag
        PL.row_start[i + 1] = (uint32_t)nnz;
    }

    PL.nnz = nnz;
    PL.col.resize(nnz);

    size_t pos = 0;
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j : rowL[i]) {
            PL.col[pos++] = j;
        }
        PL.col[pos++] = i; // diag last
    }
    assert(pos == nnz);
}

/* Build A = M + alpha S (same pattern)  */

void build_A_from_M_S(const CSRMatrix& M, const CSRMatrix& S,
                      double alpha, CSRMatrix& A)
{
    assert(M.rows == S.rows && M.cols == S.cols);
    assert(M.nnz == S.nnz);
    assert(M.row_start == S.row_start);
    assert(M.col == S.col);

    A.symmetric = true;
    A.rows = M.rows;
    A.cols = M.cols;
    A.nnz  = M.nnz;

    A.row_start = M.row_start;
    A.col       = M.col;

    A.data.resize(A.nnz);
    for (size_t k = 0; k < A.nnz; ++k) {
        A.data[k] = M.data[k] + alpha * S.data[k];
    }
}

/*  Numeric Cholesky : Algo 5.7 */

void cholesky_factorize_rowwise(const CSRMatrix& A,
                                const CSRPattern& PL,
                                CSRMatrix& L)
{
    const uint32_t n = (uint32_t)A.rows;
    assert(A.rows == A.cols);
    assert(PL.rows == A.rows && PL.cols == A.cols);

    L.symmetric = true;
    L.rows = A.rows;
    L.cols = A.cols;
    L.nnz  = PL.nnz;

    // Pattern pointers refer to PL's arrays (PL must outlive L)
    L.row_start = PL.row_start.data;
    L.col       = PL.col.data;

    L.data.resize(L.nnz);
    for (size_t k = 0; k < L.nnz; ++k) L.data[k] = 0.0;

    // Copy A into L where A has entries (A pattern is subset of PL)
    for (uint32_t i = 0; i < n; ++i) {
        size_t a0 = (size_t)A.row_start[i];
        size_t a1 = (size_t)A.row_start[i + 1];
        for (size_t ka = a0; ka < a1; ++ka) {
            uint32_t j = A.col[ka];
            size_t kl  = find_pos_in_row_binary(L, i, j);
            L.data[kl] = A.data[ka];
        }
    }

    // Factorize
    for (uint32_t i = 0; i < n; ++i) {
        size_t start   = (size_t)L.row_start[i];
        size_t diagpos = (size_t)L.row_start[i + 1] - 1;
        assert(L.col[diagpos] == i);

        // Off-diagonal entries
        for (size_t kp = start; kp < diagpos; ++kp) {
            uint32_t j = L.col[kp]; // j < i

            double s = dot_row_row_until_j(L, i, j);

            size_t diag_j = (size_t)L.row_start[j + 1] - 1;
            assert(L.col[diag_j] == j);
            double Ljj = L.data[diag_j];
            assert(Ljj > 0.0);

            L.data[kp] = (L.data[kp] - s) / Ljj;
        }

        // Diagonal
        double d = L.data[diagpos];
        for (size_t kp = start; kp < diagpos; ++kp) {
            double lij = L.data[kp];
            d -= lij * lij;
        }
        assert(d > 0.0);
        L.data[diagpos] = std::sqrt(d);
    }
}

/* Column adjacency for L^T solve */

void build_col_adjacency_lower(const CSRMatrix& L,
                               std::vector<std::vector<ColAdjEntry>>& col_adj)
{
    const uint32_t n = (uint32_t)L.rows;
    col_adj.clear();
    col_adj.resize(n);

    for (uint32_t r = 0; r < n; ++r) {
        size_t start   = (size_t)L.row_start[r];
        size_t diagpos = (size_t)L.row_start[r + 1] - 1;
        assert(L.col[diagpos] == r);

        for (size_t k = start; k < diagpos; ++k) {
            uint32_t c = L.col[k]; // c < r
            col_adj[c].push_back(ColAdjEntry{r, k});
        }
    }
}

/* Triangular solves */

void forward_solve_lower_csr(const CSRMatrix& L,
                             const double* b, double* y)
{
    const uint32_t n = (uint32_t)L.rows;

    for (uint32_t i = 0; i < n; ++i) {
        double s = b[i]; 
        size_t start   = (size_t)L.row_start[i];
        size_t diagpos = (size_t)L.row_start[i + 1] - 1;
        assert(L.col[diagpos] == i);

        for (size_t k = start; k < diagpos; ++k) {
            uint32_t j = L.col[k]; // j < i
            s -= L.data[k] * y[j];
        }
        y[i] = s / L.data[diagpos];
    }
}

void backward_solve_lowerT_with_coladj(const CSRMatrix& L,
                                       const std::vector<std::vector<ColAdjEntry>>& col_adj,
                                       const double* y, double* x)
{
    const int n = (int)L.rows;
    assert((int)col_adj.size() == n);

    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int ii = n - 1; ii >= 0; --ii) {
        uint32_t i = (uint32_t)ii;

        // x[i] -= sum_{r>i} L(r,i) * x[r]
        for (const auto& e : col_adj[i]) {
            x[i] -= L.data[e.k] * x[e.row];
        }

        size_t diagpos = (size_t)L.row_start[i + 1] - 1;
        assert(L.col[diagpos] == i);
        x[i] /= L.data[diagpos];
    }
}