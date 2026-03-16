#include "cholesky.h"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#include <vector>
#include <algorithm>

/*
 * ALGORITHM 4.2 — Build elimination tree with path compression
 *
 * Book uses 1-based indices; we use 0-based throughout.
 * Sentinel value for "no parent / not set" = N.
 *
 * Input : A stored lower-triangular (symmetric=true).
 *         Off-diagonal entries in row i have col[k] < i.
 * Output: parent[i] filled.
 */
void CholeskySolver::build_etree(const CSRMatrix &A)
{
    N = A.rows;
    parent.resize(N);
    TArray<uint32_t> ancestor(N);

    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        parent[i]   = (uint32_t)N;  /* N = "no parent"  */
        ancestor[i] = (uint32_t)N;  /* N = "not set"    */
    }

    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        /* Loop over sub-diagonal entries j < i in row i of A */
        for (uint32_t ptr = A.row_start[i]; ptr < A.row_start[i + 1]; ++ptr) {
            uint32_t j = A.col[ptr];
            if (j >= i) continue;          /* skip diagonal and upper part */

            /* Trace from j upward, applying path compression */
            uint32_t jroot = j;
            while (ancestor[jroot] != (uint32_t)N &&
                   ancestor[jroot] != i) {
                uint32_t l      = ancestor[jroot];
                ancestor[jroot] = i;       /* path compression */
                jroot           = l;
            }

            /* If jroot has no ancestor yet, set i as its parent */
            if (ancestor[jroot] == (uint32_t)N) {
                ancestor[jroot] = i;
                parent[jroot]   = i;
            }
        }
    }
}

/* 
 * ALGORITHM 4.3 — Compute row sparsity patterns of L
 *
 * For each row i, we walk up the elimination tree from every sub-diagonal
 * entry k in row i of A, collecting ancestors until we hit one already
 * marked for row i.
 *
 * Result is stored in L_pat (owns the memory).
 * Columns within each row are sorted ascending; diagonal is last.
 * L.row_start and L.col are set to point into L_pat afterwards.
 */
void CholeskySolver::build_pattern(const CSRMatrix &A)
{
    /* mark[j] == i  ⟺  column j already added to row i */
    TArray<uint32_t> mark(N);
    for (uint32_t i = 0; i < (uint32_t)N; ++i) mark[i] = (uint32_t)N;

    /* Collect column lists per row */
    std::vector<std::vector<uint32_t>> rows(N);

    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        /* Diagonal is always in L */
        mark[i] = i;
        rows[i].push_back(i);

        /* Sub-diagonal entries of A in row i */
        for (uint32_t ptr = A.row_start[i]; ptr < A.row_start[i + 1]; ++ptr) {
            uint32_t k = A.col[ptr];
            if (k >= i) continue;

            /* Walk up the elimination tree from k */
            uint32_t j = k;
            while (j < (uint32_t)N && mark[j] != i) {
                mark[j] = i;             /* flag as seen for row i */
                rows[i].push_back(j);    /* add to pattern         */
                j = parent[j];           /* move up the tree       */
            }
        }

        /* Sort ascending (diagonal i is max, ends up last) */
        std::sort(rows[i].begin(), rows[i].end());
    }

    /* Build L_pat from rows[] */
    L_pat.symmetric = false;
    L_pat.rows = N;
    L_pat.cols = N;
    L_pat.row_start.resize(N + 1);
    L_pat.row_start[0] = 0;
    for (uint32_t i = 0; i < (uint32_t)N; ++i)
        L_pat.row_start[i + 1] = L_pat.row_start[i] + (uint32_t)rows[i].size();

    size_t nnz_L = L_pat.row_start[N];
    L_pat.nnz = nnz_L;
    L_pat.col.resize(nnz_L);
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        uint32_t base = L_pat.row_start[i];
        for (uint32_t k = 0; k < (uint32_t)rows[i].size(); ++k)
            L_pat.col[base + k] = rows[i][k];
    }

    /* Wire up L's raw pointers to L_pat's owned arrays */
    L.symmetric = false;
    L.rows      = N;
    L.cols      = N;
    L.nnz       = nnz_L;
    L.row_start = L_pat.row_start.data;
    L.col       = L_pat.col.data;
    L.data.resize(nnz_L);

    printf("[Cholesky] N=%zu  nnz(A)=%u  nnz(L)=%zu  fill=%.2fx\n",
           N, A.row_start[N], nnz_L,
           (double)nnz_L / (double)A.row_start[N]);

    /* Verify each row is sorted ascending with diagonal last */
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        for (uint32_t p = L_pat.row_start[i]; p+1 < L_pat.row_start[i+1]; ++p) {
            if (L_pat.col[p] >= L_pat.col[p+1]) {
                printf("[Cholesky] ERROR: row %u not sorted: col[%u]=%u >= col[%u]=%u\n",
                       i, p, L_pat.col[p], p+1, L_pat.col[p+1]);
                assert(false);
            }
        }
        if (L_pat.col[L_pat.row_start[i+1]-1] != i) {
            printf("[Cholesky] ERROR: row %u diagonal not last\n", i);
            assert(false);
        }
    }
    printf("[Cholesky] Pattern OK.\n"); fflush(stdout);
}

/* 
 * ALGORITHM 5.6 — Sparse right-looking Cholesky
 *
 * For j = 0 .. N-1  (process column j):
 *   1.  l_jj = sqrt(l_jj)
 *   2.  l_ij /= l_jj   for all i > j in col(L, j)
 *   3.  For k > j in col(L, j):
 *         For i >= k   in col(L, j):
 *           l_ik -= l_ij * l_kj
 *
 * col(L, j) = set of rows that have a nonzero in column j.
 * In row-CSR storage we precompute this as col_rows[j].
 *
 * find_ptr(i, j) returns the index into L.data for entry (i,j).
 * Columns within each row are sorted, so binary search is used.
 */
void CholeskySolver::numeric(const CSRMatrix &A)
{
    /* Sanity checks */
    printf("[Cholesky] numeric: N=%zu nnz=%zu row_start[N]=%u\n",
           N, L.nnz, L.row_start[N]);
    fflush(stdout);
    assert(L.row_start[N] == (uint32_t)L.nnz);
    assert(L.data.size == L.nnz);

    /* Zero-initialise L.data */
    memset(L.data.data, 0, L.nnz * sizeof(double));

    /* Binary search: find data index for (row i, col j) in L */
    auto find_ptr = [&](uint32_t i, uint32_t j) -> uint32_t {
        uint32_t lo = L.row_start[i];
        uint32_t hi = L.row_start[i + 1];
        while (lo < hi) {
            uint32_t mid = (lo + hi) >> 1;
            if      (L.col[mid] == j) return mid;
            else if (L.col[mid] <  j) lo = mid + 1;
            else                      hi = mid;
        }
        printf("[Cholesky] ERROR: entry (%u,%u) not found in L pattern!\n", i, j);
        fflush(stdout);
        /* print row i's columns for debugging */
        printf("  Row %u columns:", i);
        for (uint32_t p = L.row_start[i]; p < L.row_start[i+1]; ++p)
            printf(" %u", L.col[p]);
        printf("\n");
        fflush(stdout);
        assert(false);
        return 0;
    };

    printf("[Cholesky] Copying A into L...\n"); fflush(stdout);
    /* Copy lower triangle of A into L */
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        for (uint32_t ptr = A.row_start[i]; ptr < A.row_start[i + 1]; ++ptr) {
            uint32_t j = A.col[ptr];
            if (j > i) continue;
            /* Verify (i,j) is in L pattern before calling find_ptr */
            uint32_t lo = L.row_start[i];
            uint32_t hi = L.row_start[i + 1];
            bool found = false;
            for (uint32_t p = lo; p < hi; ++p) {
                if (L.col[p] == j) { found = true; break; }
            }
            if (!found) {
                printf("[Cholesky] MISSING in L: A(%u,%u) not in L pattern!\n", i, j);
                fflush(stdout);
                /* print L row i */
                printf("  L row %u:", i);
                for (uint32_t p = lo; p < hi; ++p) printf(" %u", L.col[p]);
                printf("\n"); fflush(stdout);
                assert(false);
            }
            L.data.data[find_ptr(i, j)] = A.data.data[ptr];
        }
    }

    printf("[Cholesky] Building col_rows...\n"); fflush(stdout);
    /* Precompute col_rows[j] = sorted list of row indices i > j */
    std::vector<std::vector<uint32_t>> col_rows(N);
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        for (uint32_t ptr = L.row_start[i]; ptr < L.row_start[i + 1]; ++ptr) {
            uint32_t j = L.col[ptr];
            if (j < i) col_rows[j].push_back(i);
        }
    }

    printf("[Cholesky] Starting factorization loop...\n"); fflush(stdout);
    for (uint32_t j = 0; j < (uint32_t)N; ++j) {

        /* Step 1: l_jj = sqrt(l_jj) */
        uint32_t diag_ptr = find_ptr(j, j);
        double   ljj      = L.data.data[diag_ptr];
        if (ljj <= 0.0) {
            printf("[Cholesky] Non-positive pivot at j=%u: ljj=%.6e\n", j, ljj);
            assert(ljj > 0.0 && "Non-positive pivot: matrix is not SPD");
        }
        ljj = sqrt(ljj);
        L.data.data[diag_ptr] = ljj;

        /* Step 2: scale sub-diagonal entries in column j */
        for (uint32_t i : col_rows[j]) {
            L.data.data[find_ptr(i, j)] /= ljj;
        }

        /* Step 3: rank-1 update of trailing sub-matrix */
        const auto &cj = col_rows[j];          /* rows with nonzero in col j */
        for (size_t ki = 0; ki < cj.size(); ++ki) {
            uint32_t k   = cj[ki];
            double   lkj = L.data.data[find_ptr(k, j)];

            /* update diagonal l_kk */
            L.data.data[find_ptr(k, k)] -= lkj * lkj;

            /* update l_ik for i > k also in col_rows[j] */
            for (size_t ii = ki + 1; ii < cj.size(); ++ii) {
                uint32_t i   = cj[ii];
                double   lij = L.data.data[find_ptr(i, j)];
                L.data.data[find_ptr(i, k)] -= lij * lkj;
            }
        }
    }
}


void CholeskySolver::factorize(const CSRMatrix &A)
{
    printf("[Cholesky] Starting build_etree...\n"); fflush(stdout);
    build_etree(A);
    printf("[Cholesky] Starting build_pattern...\n"); fflush(stdout);
    build_pattern(A);
    printf("[Cholesky] Starting numeric...\n"); fflush(stdout);
    numeric(A);
    printf("[Cholesky] Done.\n"); fflush(stdout);
}

/* 
 * Forward substitution: solve  L y = b
 * Row i:  y[i] = ( b[i] - sum_{j < i} L(i,j) * y[j] ) / L(i,i)
 */
void CholeskySolver::forward_solve(const double *b, double *y) const
{
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        double s   = b[i];
        double lii = 1.0;
        for (uint32_t ptr = L.row_start[i]; ptr < L.row_start[i + 1]; ++ptr) {
            uint32_t j = L.col[ptr];
            if      (j < i) s   -= L.data.data[ptr] * y[j];
            else if (j == i) lii = L.data.data[ptr];
        }
        y[i] = s / lii;
    }
}

/* 
 * Backward substitution: solve  L^T x = y
 * Process rows bottom-to-top.
 * For i = N-1 .. 0:
 *   x[i] /= L(i,i)
 *   for j < i with L(i,j) != 0:  x[j] -= L(i,j) * x[i]
*/
void CholeskySolver::backward_solve(const double *y, double *x) const
{
    memcpy(x, y, N * sizeof(double));

    for (int i = (int)N - 1; i >= 0; --i) {
        double lii = 1.0;
        for (uint32_t ptr = L.row_start[i]; ptr < L.row_start[i + 1]; ++ptr)
            if (L.col[ptr] == (uint32_t)i) { lii = L.data.data[ptr]; break; }

        x[i] /= lii;

        /* Scatter x[i] into x[j] for j < i (transpose of L row i) */
        for (uint32_t ptr = L.row_start[i]; ptr < L.row_start[i + 1]; ++ptr) {
            uint32_t j = L.col[ptr];
            if (j < (uint32_t)i)
                x[j] -= L.data.data[ptr] * x[i];
        }
    }
}

/* A x = b */
void CholeskySolver::solve(const double *b, double *x) const
{
    TArray<double> y(N);
    forward_solve(b, y.data);
    backward_solve(y.data, x);
}