#include <stdio.h>
#include "mesh.h"
#include "sphere.h"
#include "P1.h"
#include "sparse_matrix.h"
#include "sparse_cholesky.h"

int main(){
	CSRPattern pattern {0};
	// pattern from fig 4.2 in the book
	pattern.rows = 8;
	pattern.cols = 8;
	pattern.nnz = 24;
	uint32_t row_start[] = {0,3,7,10,13,16,18,20,24};
	uint32_t col[] = {
		0,4,5,
		1,3,4,7,
		2,3,7,
		1,2,3,
		0,1,4,
		0,5,
		6,7,
		1,2,6,7};
	pattern.row_start.data = row_start;
	pattern.row_start.size = 9;
	pattern.row_start.capacity = 9;

	pattern.col.data = col;
	pattern.col.size = 24;
	pattern.col.capacity = 24;
	// Mesh m;
	// load_sphere(m,0);
	// build_P1_CSRPattern(m,pattern);
	CSRMatrix A(pattern,0);
	TArray<int32_t> etree;
	elimination_tree(A,etree);
	A.print_pattern();
	printf("\n");
	for (uint32_t i = 0; i<etree.size; i++){
		printf("%u : %d\n",i+1,etree[i]+1);
	}
}
