#include "sparse_tri_solve.h"

/*
 * Solve a lower triangular system were b and x are full vector
 */
void sparse_tri_solve(CSRMatrix &L, double *__restrict b,double *__restrict x){
	for (size_t row = 0; row < L.rows; row++){
		double sum = 0;
		for(uint32_t k = L.row_start[row]; k<L.row_start[row+1];k++){
			uint32_t col = L.col[k];
			if (col == row){
				x[row] = (b[row] - sum)/L.data[k];
				break;
			}
			sum += L.data[k] * x[col];
		}
	}
}

/*
 * Solve a lower triangular system were b and x are sparse vector
 * Assumes L has non-zero coeficients on its diagonal
 */
void sparse_tri_solve(CSRMatrix &L,
		double *b, uint32_t *b_ind, size_t b_size,
		double *x, uint32_t *x_ind, size_t x_size,
		uint32_t i){
	size_t current_b_k = 0;
	for (size_t row_i = 0; row_i < x_size; row_i++){
		size_t row = x_ind[row_i];
		if (row >= i) return ;
		double sum = 0;
		size_t current_x_k = 0;
		for(uint32_t k = L.row_start[row]; k<L.row_start[row+1];k++){
			uint32_t col = L.col[k];
			if (col == row){
				while (b_ind[current_b_k]<row && current_b_k < b_size-1) current_b_k++;
				if (b_ind[current_b_k] == row){
					x[row_i] = (b[current_b_k] - sum)/L.data[k];
				}
				else {
					x[row_i] = - sum/L.data[k];
				}
				break;
			}
			while (x_ind[current_x_k]<col && current_x_k < x_size-1) current_x_k++;
			if (x_ind[current_x_k] == col){
				sum += L.data[k] * x[current_x_k];
			}
		}
	}
}
