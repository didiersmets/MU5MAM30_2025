#include <stdio.h>
#include <string.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness.h"
#include "algorithm"

/* CSRMatrix variants */

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	/* Your implementation goes here.
	 * Use a VTAdjacency structure (see include/matrix/adjacency.h)
	 */

	 /*
		build the sparse matrix pattern to store information on how 
		all the vertixes are connected.
		The main idea is to use the Adjency pattern built before:
		
		for the row array we can just use the offset array we computed
		in the constructor of Adjency since it already tells us on 
		each row how many entries we have, since in each row we will have
		as many entries as the degree of that specific vtx.

		for the column array we instead need to loop over the vtri array
		in this array we will have informations on which are the vtxs to which
		our vtx is connected to, for each vtx we need to loop from
		vtri[ offset[vtx] ] to vtri[ offset[vtx] + degree ]

	 */
	P.symmetric = true;
	P.rows = m.vertex_count();
	P.cols = m.vertex_count();
	VTAdjacency adj(m);

	P.row_start.resize( P.rows + 1 );
	P.col.resize( P.cols + 3 * m.triangle_count() ); 
	
	int nnz = 0;
	for( int a = 0; a < P.rows; a++)
	{
		//also try with offset array in adjacency
		// P.row_start[a] = adj.offset[a]
		P.row_start[a] = nnz;

		int init_nnz = nnz;

		int start = adj.offset[a];
		int stop = start + adj.degree[a];

		for( int i = start; i < stop; i++ )
		{
			//retrieve the entries from the vtri
			int b = adj.vtri[i].next;
			int c = adj.vtri[i].prev;
			
			/*
			since we are iterating for the number of a connected to vtx
			and not on the triangles that vtx is part of we will end up with 
			redundant entries in the column array, to avoid that we need to see
			if the entries are already present 
			*/ 
			bool b_present = false;
			bool c_present = false;
			int tmp_add_elem = nnz - init_nnz;
			
			for( int j = 0; j < tmp_add_elem; j++ )
			{
				if( P.col[init_nnz + j] == b )
					b_present = true;
				if( P.col[init_nnz + j] == c )
					c_present = true;
			}
			// since the mtx is symmetric we only add the lower diagonal elements
			// so we check b < a and c < a
			if( b < a && !b_present )
				P.col[nnz++] = b;
			if( c < a && !c_present )
				P.col[nnz++] = c;
			
		}
		P.col[nnz++] = a; // add the diagonal entry
	}

	P.col.resize(nnz);
	P.col.shrink_to_fit();

	// reorder the column indices in each row in increasing order
	for (size_t a = 0; a < P.rows; ++a) {
		uint32_t *start = &P.col[P.row_start[a]];
		uint32_t *end = &P.col[P.row_start[a + 1]];
		std::sort(start, end);
	}

	 
}

void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	assert(P.row_start.size == vtx_count + 1);

	M.symmetric = true;
	M.rows = M.cols = vtx_count;
	M.nnz = P.col.size;
	M.row_start = P.row_start.data;
	M.col = P.col.data;
	M.data.resize(M.nnz);
	for (size_t i = 0; i < M.nnz; ++i) {
		M.data[i] = 0.0;
	}

	/* Your implementation goes here */

	/*
	each row and col if the matrix represent a specific basis function
	we need to loop over all the triangles and compute the local mass matrix
	for each triangle we then need to add the contributions to the global mass matrix
	*/
	const TArray<uint32_t> &idx = m.indices;
	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = { (double)B[0] - (double)A[0], (double)B[1] - (double)A[1], (double)B[2] - (double)A[2] };
		Vec3d AC = { (double)C[0] - (double)A[0], (double)C[1] - (double)A[1], (double)C[2] - (double)A[2] };
		// as discussed in mass.h we only need to store 2 values for the local mass matrix
		double Mloc[2]; 
		//create the local mass matrix
		mass(AB, AC, Mloc);
		double diagonal_term = Mloc[0];
		double off_diagonal_term = Mloc[1];
		// add contributions to global mass matrix
		M(a, a) += diagonal_term;
		M(b, b) += diagonal_term;
		M(c, c) += diagonal_term;
		// off diagonal terms
		const unsigned int ab_max = std::max(a, b);
		const unsigned int ab_min = std::min(a, b);
		M(ab_max, ab_min) += off_diagonal_term;

		const unsigned int bc_max = std::max(b, c);
		const unsigned int bc_min = std::min(b, c);
		M(bc_max, bc_min) += off_diagonal_term;

		const unsigned int ca_max = std::max(c, a);
		const unsigned int ca_min = std::min(c, a);
		M(ca_max, ca_min) += off_diagonal_term;
	}
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	assert(P.row_start.size == vtx_count + 1);

	S.symmetric = true;
	S.rows = S.cols = vtx_count;
	S.nnz = P.col.size;
	S.row_start = P.row_start.data;
	S.col = P.col.data;
	S.data.resize(S.nnz);
	for (size_t i = 0; i < S.nnz; ++i) {
		S.data[i] = 0.0;
	}

	/* Your implementation goes here */
	/*
	same logic as per the stiffness matrix but now we compute the stiffness matrix
	*/
	const TArray<uint32_t> &idx = m.indices;
	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = { (double)B[0] - (double)A[0], (double)B[1] - (double)A[1], (double)B[2] - (double)A[2] };
		Vec3d AC = { (double)C[0] - (double)A[0], (double)C[1] - (double)A[1], (double)C[2] - (double)A[2] };
		// as discussed in stiffness.h we only need to store 2 values for the local stiffness matrix
		double Sloc[6]; 
		//create the local stiffness matrix
		stiffness(AB, AC, Sloc);
		// add contributions to global stiffness matrix
		S(a, a) += Sloc[0];
		S(b, b) += Sloc[1];
		S(c, c) += Sloc[2];
		// off diagonal terms
		const unsigned int ab_max = std::max(a, b);
		const unsigned int ab_min = std::min(a, b);
		S(ab_max, ab_min) += Sloc[3];

		const unsigned int bc_max = std::max(b, c);
		const unsigned int bc_min = std::min(b, c);
		S(bc_max, bc_min) += Sloc[4];

		const unsigned int ca_max = std::max(c, a);
		const unsigned int ca_min = std::min(c, a);
		S(ca_max, ca_min) += Sloc[5];
	}
}

/* FEMatrix variants */

void build_P1_mass_matrix(const Mesh &m, FEMatrix &M)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();

	M.fem_type = FEMatrix::P1_cst;
	M.m = &m;
	M.rows = M.cols = vtx_count;

	M.diag.resize(vtx_count);
	memset(M.diag.data, 0, vtx_count * sizeof(double));

	M.off_diag.resize(tri_count);
	const TArray<uint32_t> &idx = m.indices;
	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = { (double)B[0] - (double)A[0],
			     (double)B[1] - (double)A[1],
			     (double)B[2] - (double)A[2] };
		Vec3d AC = { (double)C[0] - (double)A[0],
			     (double)C[1] - (double)A[1],
			     (double)C[2] - (double)A[2] };
		double Mloc[2];
		mass(AB, AC, Mloc);
		M.diag[a] += Mloc[0];
		M.diag[b] += Mloc[0];
		M.diag[c] += Mloc[0];
		M.off_diag[t] = Mloc[1];
	}
}

void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();

	S.fem_type = FEMatrix::P1_sym;
	S.m = &m;
	S.rows = S.cols = vtx_count;

	S.diag.resize(vtx_count);
	memset(S.diag.data, 0, vtx_count * sizeof(double));

	S.off_diag.resize(3 * tri_count);
	const TArray<uint32_t> &idx = m.indices;
	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = { (double)B[0] - (double)A[0],
			     (double)B[1] - (double)A[1],
			     (double)B[2] - (double)A[2] };
		Vec3d AC = { (double)C[0] - (double)A[0],
			     (double)C[1] - (double)A[1],
			     (double)C[2] - (double)A[2] };
		double Sloc[6];
		stiffness(AB, AC, Sloc);
		S.diag[a] += Sloc[0];
		S.diag[b] += Sloc[1];
		S.diag[c] += Sloc[2];
		S.off_diag[3 * t + 0] = Sloc[3];
		S.off_diag[3 * t + 1] = Sloc[4];
		S.off_diag[3 * t + 2] = Sloc[5];
	}
}
