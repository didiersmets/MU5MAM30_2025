#include <math.h>
#include <algorithm>
#include "sparse_matrix.h"
#include "array.h"

void cholesky_factorization(CSRMatrix &A, CSRPattern &P, CSRMatrix &L){
    //Check that A is SPD and that the sparsity patten is the same size as A.
    assert(A.symmetric);
    assert(A.cols == P.cols);
    assert(A.rows == P.rows);

    //Initialize L
    size_t n = P.cols;
    L.cols = n;
    L.rows = n;
    
    L.data.resize(P.nnz);
    L.row_start = P.row_start.data;
    L.col = P.col.data;
    L.symmetric = false;

    //set all entries in L to 0;
    for (size_t i = 0; i < P.nnz; i++) {
		L.data[i] = 0.0;
	}

    //Fille entries of L that are alreadu filled in A
    for(size_t i=0; i <n; i ++){
        for(size_t j = A.row_start[i]; j < A.row_start[i+1]; j++){
            L(i, A.col[j]) = A.data[j];
        }
    }

    for(size_t j = 0; j < n; j++){
        for(size_t k = P.row_start[j]; k < P.row_start[j+1]; k++){ //look at the jth row
            size_t kcol = P.col[k];
            if(kcol >= j) break; //We assume that the columns are ordered
            for(size_t i = j; i < n; i++){ //now we look at the kth column, indicies bigger than j
                for(size_t l = P.row_start[i]; l < P.row_start[i+1]; l++){ 
                    //this is a bit slow but we go by row to see if (i, kcol) is in the pattern
                    if(P.col[l] == kcol){//it is! so we can add it
                        L(i, j) = L(i, j) - L.data[l]*L.data[k];
                    }
                }
            }
        }
        L.data[L.row_start[j+1]-1] = sqrt(L.data[L.row_start[j+1]-1]);
        for(size_t i = j+1; i <n; i ++){
            for(size_t l = P.row_start[i]; l < P.row_start[i+1]; l++){
                    if(P.col[l] == j){
                        L.data[l] =  L.data[l]/L.data[L.row_start[j+1]-1];
                    }
                }
            }
        }
    
}

void solve_cholesky(CSRMatrix &L,  const double *__restrict b, double *__restrict x){
    size_t n = L.cols;

    //Forward substitution
    for(size_t i = 0; i < n; i++){
        x[i] = b[i];
        for(size_t j = L.row_start[i]; j < L.row_start[i+1]-1; j++){
             x[i] -= L.data[j]*x[L.col[j]];
        }
        x[i] /= L.data[L.row_start[i+1]-1];
    }
    
    //Backward subsitution
    for (size_t i = n; i-- > 0; ) {

        double xi = x[i] / L.data[L.row_start[i+1]-1]; //store x_i as to not modify

        for (size_t k = L.row_start[i]; k < L.row_start[i+1]-1; k++) {
            size_t j = L.col[k];   // j < i
            x[j] -= L.data[k] * xi;
        }

        x[i] = xi;
    }


}
