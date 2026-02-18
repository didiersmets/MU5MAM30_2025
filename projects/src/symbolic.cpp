#include "symbolic.h"
#include "sparse_matrix.h"
#include "array.h"

//Define this function here again to be able to switch between my P1 and the solution P1
static bool find_2(const TArray<uint32_t> & arr, uint32_t target,  size_t start, size_t stop){
	assert(start <= stop && stop <= arr.size);

	for(size_t i = start; i < stop; i ++){
		if(arr[i] == target){return true;}
	}
	return false;
}

void construct_etree(CSRPattern &P, TArray<uint32_t> &parent){
    assert(P.cols == P.rows); //check that matrix is square
    assert(P.symmetric); //check that input is symmetric
    
    parent.resize(P.cols);

    size_t n = P.cols;
    TArray<uint32_t> anc(n);
    for(size_t i = 0; i<n; i++){
        // don't use 0 for initialization since indexing starts from 0
       parent[i] = n+1; 
       anc[i] = n+1;
       for(size_t j = P.row_start[i]; j < P.row_start[i+1]; j ++ ){ //loop over each row
            if(P.col[j] >= i) continue; //not necessary per say but just to be safe
            uint32_t jroot = P.col[j];
            while(anc[jroot] != n+1 && anc[jroot] != i){
                uint32_t l = anc[jroot];
                anc[jroot] = i; // path compression step
                jroot = l;
            }
            if(anc[jroot] == n+1){
                anc[jroot] = i;
                parent[jroot] = i;
            }
       }
    }
}

void construct_L_sparsity_pattern(CSRPattern &A, CSRPattern &L, TArray<uint32_t> &parent){
    assert(A.cols == A.rows); //check that matrix is square
    assert(A.symmetric); //check that input is symmetric
    assert(A.cols <= parent.size); //check that we will fit the data in the parent vector
    
    size_t n = A.cols;
    L.rows = n;
    L.cols = n;
    L.symmetric = false;

    L.row_start.resize(n +1);
    L.row_start[0] = 0;
    L.col.resize(n*n/2 + n); // initalize to max number of nnz

    TArray<uint32_t> mark(n);
    size_t nnz = 0;

    for(size_t i = 0; i < n ; i ++){
        mark[i] = i;
        for(size_t k = A.row_start[i]; k < A.row_start[i+1]; k ++){
            size_t j = A.col[k];
            while(mark[j] != i){
                mark[j] = i;
                //if we did not already add j, we add it
                if(!find_2(L.col, j, L.row_start[i], nnz)){
                    L.col[nnz++] = j;
                }
                j = parent[j];
                if (j >= n) break;
            }
            
        }
        L.col[nnz++] = i;
        L.row_start[i+1] = nnz;
    }

    L.nnz = nnz;
    L.col.resize(L.nnz);
    L.col.shrink_to_fit();

    //reorder the indicies
	uint32_t stop = 0;
	for(size_t i = 0; i < n; i ++){
		uint32_t start = stop;
		stop = L.row_start[i+1];

		for(size_t j = start+1; j < stop; j++){
			for(size_t k = start; k <j; k++){
				if(L.col[j] < L.col[k]){
					uint32_t tmp = L.col[j];
					L.col[j] = L.col[k];
					L.col[k] = tmp;
				}
			}
		}
	}

}
