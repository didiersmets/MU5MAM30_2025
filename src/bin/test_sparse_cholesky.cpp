#include <stdio.h>
#include "mesh.h"
#include "sphere.h"
#include "P1.h"
#include "sparse_matrix.h"
#include "sparse_cholesky.h"
#include "sparse_tri_solve.h"

int main(){
	CSRPattern pattern {0};
# if 1
	CSRMatrix A;
	A.rows = 5;
	A.cols = 5;
	A.nnz = 24;
	uint32_t row_start[] = {0,1,5,7,10,13};
	uint32_t col[] = {
		0,
		1,2,3,4,
		1,2,
		1,3,4,
		1,3,4,
		};
	double data[] = {4,2,1,1,1,1,2,1,2,1,1,1,2};
	A.row_start = row_start;
	A.col = col;
	A.data.data = data;
	A.data.size = sizeof(data)/sizeof(data[0]);
	A.data.capacity = sizeof(data)/sizeof(data[0]);
#endif
# if 0
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
	CSRMatrix A(pattern,1);
#endif
# if 0
	Mesh m;
	load_sphere(m,0);
	build_P1_CSRPattern(m,pattern);
	CSRMatrix A(pattern,1);
	build_P1_mass_matrix(m,A);
#endif
	A.print();
	double b[8] = {1,2,3,4,5,6,7,8};
	uint32_t b_ind[8] = {0,1,2,3,4,5,6,7};
	double x[8] = {0};

	TArray<uint32_t> etree;
	elimination_tree(A,etree);
	A.print_pattern();
	printf("\n");
	for (uint32_t i = 0; i<etree.size; i++){
		printf("%u : %u\n",i,etree[i]);
	}
	CSRPattern cho_patt {0};
	cholesky_sparsity_pattern(A,etree,cho_patt);
	CSRMatrix L(cho_patt,1);
	L.print_pattern();
	sparse_tri_solve(L,b,b_ind,8,x,b_ind,8,2);
	for (uint32_t i = 0; i<etree.size; i++){
		printf("%u : %lf\n",i,x[i]);
	}

	CSRMatrix L2(cho_patt,0);
	sparse_cholesky(A,L2);
	L2.print();
}
