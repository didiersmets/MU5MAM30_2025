#include "sparse_cholesky.h"

void myqsort(TArray<uint32_t> &a, size_t start, size_t end){
	if (start+1 >= end) return;

	uint32_t pivot = a[start];
	size_t mid = start+1;
	for (size_t i=start+1;i<end;i++){
		if (a[i]<=pivot){
			uint32_t temp = a[i];
			a[i] = a[mid];
			a[mid] = temp;
			mid++;
		}
	}
	a[start] = a[mid-1];
	a[mid-1] = pivot;

	myqsort(a,start,mid-1);
	myqsort(a,mid,end);
}

void elimination_tree(const CSRMatrix &A,TArray<uint32_t> &parent){
	parent.resize(A.rows);
	TArray<uint32_t> ancestor(A.rows,UINT32_MAX);
	for (uint32_t i = 0; i<A.rows; i++){
		parent[i] = UINT32_MAX;
		for(uint32_t k = A.row_start[i]; k<A.row_start[i+1];k++){
			uint32_t j = A.col[k];
			if (j >= i) break;
			while (ancestor[j] != UINT32_MAX && ancestor[j] != i){
				uint32_t l = ancestor[j];
				ancestor[j] = i;
				j = l;
			}
			if (ancestor[j] == UINT32_MAX){
				ancestor[j] = i;
				parent[j] = i;
			}
		}

	}
}

void cholesky_sparsity_pattern(const CSRMatrix &A,const TArray<uint32_t> &parent,CSRPattern &cho_patt){
	cho_patt.rows = A.rows;
	cho_patt.cols = A.cols;
	cho_patt.row_start.resize(0);
	cho_patt.row_start.push_back(0);
	cho_patt.col.resize(0);

	TArray<TArray<uint32_t>> rowL(A.rows);
	for (uint32_t i = 0; i< A.rows; i++){
		rowL[i].size = 0;
		rowL[i].capacity = 0;
		rowL[i].data = nullptr;
	}

	TArray<uint32_t> mark(A.rows,UINT32_MAX);
	for (uint32_t i = 0; i < A.rows; i++){
		rowL[i].push_back(i);
		mark[i] = i;
		for(uint32_t k = A.row_start[i]; k<A.row_start[i+1];k++){
			uint32_t j = A.col[k];
			if (j >= i) break;
			while (mark[j] != i){
				mark[j] = i;
				rowL[i].push_back(j);
				j = parent[j];
			}
		}
	}
	for (uint32_t i = 0; i<A.rows ; i++){
		cho_patt.row_start.push_back(rowL[i].size + cho_patt.row_start[i]);
	}
	cho_patt.nnz = cho_patt.row_start[A.rows];
	cho_patt.col.reserve(cho_patt.nnz);
	for (uint32_t i = 0; i<A.rows ; i++){
		for (uint32_t j=0;j<rowL[i].size;j++){
			cho_patt.col.push_back(rowL[i][j]);
		}
	}
	for(size_t i = 0;i<cho_patt.rows;i++) myqsort(cho_patt.col,cho_patt.row_start[i],cho_patt.row_start[i+1]);
}
void sparse_cholesky(CSRMatrix &A,CSRMatrix &L){
	// TArray<uint32_t> etree;
	// elimination_tree(A,etree);
	// CSRPattern cho_patt;
	// cholesky_sparsity_pattern(A,etree,cho_patt);
	// CSRMatrix L(cho_patt,0);
	L(0,0) = sqrt(A(0,0));
	for (uint32_t i = 1; i<A.rows ;i++){
		uint32_t b_start = A.row_start[i];
		double *b = A.data.data + b_start;
		uint32_t *b_ind = A.col + b_start;
		size_t b_size = 0;
		while (b_ind[b_size]<i && b_start + b_size < A.row_start[i+1]) b_size++;
		uint32_t x_start = L.row_start[i];
		double *x = L.data.data + x_start;
		uint32_t *x_ind = L.col + x_start;
		size_t x_size = 0;
		while (x_ind[x_size]<i && x_start + x_size < L.row_start[i+1]) x_size++;
		sparse_tri_solve(L,b,b_ind,b_size,x,x_ind,x_size,i);
		double dot_product = blas_dot(x,x,x_size);
		L.data[x_start+x_size] = sqrt(A.data[b_start+b_size] - dot_product);
	}
}
