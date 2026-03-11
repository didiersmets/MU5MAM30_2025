#include "sparse_cholesky.h"

void elimination_tree(const CSRMatrix &A,TArray<int32_t> &parent){
	parent.resize(A.rows);
	TArray<int32_t> ancestor(A.rows,-1);
	for (uint32_t i = 0; i<A.rows; i++){
		parent[i] = -1;
		for(uint32_t k = A.row_start[i]; k<A.row_start[i+1];k++){
			uint32_t j = A.col[k];
			if (j >= i) break;
			while (ancestor[j] != -1 && ancestor[j] != (int32_t) i){
				uint32_t l = ancestor[j];
				ancestor[j] = i;
				j = l;
			}
			if (ancestor[j] == -1){
				ancestor[j] = i;
				parent[j] = i;
			}
		}

	}
}
