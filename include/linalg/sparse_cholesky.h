#include "array.h"
#include "sparse_matrix.h"
#include "math.h"
#include "sparse_tri_solve.h"
#include "tiny_blas.h"

void elimination_tree(const CSRMatrix &A,TArray<uint32_t> &parent);
void cholesky_sparsity_pattern(const CSRMatrix &A,const TArray<uint32_t> &parent,CSRPattern &cho_patt,CSRPattern &cho_anti_diag_patt);
void sparse_cholesky(CSRMatrix &A,CSRMatrix &L);
void anti_transpose(CSRMatrix &L1, CSRMatrix &L2);

struct CholeskySolver {
	TArray<uint32_t> etree;
	CSRPattern cho_patt;
	CSRPattern cho_anti_diag_patt;
	CSRMatrix L;
	CSRMatrix L_anti_transpose;
	TArray<double> tmp_buffer;

	CholeskySolver();
	CholeskySolver(CSRMatrix &M);
	void update_same_pattern(CSRMatrix &M);
	void solve(double *b,double *x);
};
