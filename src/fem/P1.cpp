#include <stdio.h>
#include <string.h>
#include <algorithm>

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
#include "math_utils.h"

/* CSRMatrix variants */

// Auxiliary function
bool is_in_first_idx(uint32_t target, uint32_t *array, size_t nb_idx_check)
{
	for (size_t i=0; i<nb_idx_check; i++) {
		if (array[i]==target) {
			return true;
		}
	}
	return false;
}

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	/* Your implementation goes here.
	 * Use a VTAdjacency structure (see include/mesh/adjacency.h)
	 */
	VTAdjacency vtadj(m);
	P.symmetric = true;
	P.rows = m.vertex_count();
	P.cols = m.vertex_count();


	/* Note that each vertex index i correspond to a line in the sparse matrix 
	   therefore we build the pattern line by line */
	P.col.resize(3*m.triangle_count() + P.rows); // 3 sides per triangle + diag term

	size_t nnz = 0;
	for (size_t i = 0; i<P.rows; i++) {
		P.row_start[i] = nnz; // The i line begins at this nnz-th value
		size_t line_nnz = 0;
		uint32_t *line_start = &P.col[nnz];

		/* iterate on every triangle connected to vertex i */
		uint32_t tri_start = vtadj.offset[i]; // first triangle in VTAdjency structure
		uint32_t tri_stop = tri_start + vtadj.degree[i]; // next to last triangle

		for (size_t tri_index = tri_start; tri_index<tri_stop; tri_index++) { // triangle ijk
			uint32_t j = vtadj.vtri[tri_index].next;
			uint32_t k = vtadj.vtri[tri_index].prev;

			if (j<i && !is_in_first_idx(j, line_start, line_nnz)) { // Check if j was not already encounter
				nnz++;
				line_nnz++;
				P.col[nnz] = j;
			}

			if (k<i && !is_in_first_idx(k, line_start, line_nnz)) { // Check if k was not already encounter
				nnz++;
				line_nnz++;
				P.col[nnz] = k;
			}
		}
		nnz++;
		P.col[nnz] = i; // i is always connected to himself (diag term)
	}
	P.row_start[P.rows] = nnz;
	P.nnz = nnz;
	P.col.resize(nnz); // last size was upper estimate

	/* each col needs to be sorted in the final pattern*/
	for (size_t i = 0; i < P.rows; i++) {
		size_t line_start = P.row_start[i];
		size_t line_end = P.row_start[i+1];

		std::sort(&P.col[line_start], &P.col[line_end]);
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

	/* Since we can compute the coefficients for each triangle, we will build M by
	by additionning the contribution of each triangle */

	for (size_t tri_index = 0; tri_index<tri_count; tri_index+=3) { // Triangle ABC d'indices i,j,k
		uint32_t i = m.indices[tri_index];
		uint32_t j = m.indices[tri_index+1];
		uint32_t k = m.indices[tri_index+2];

		Vec3 A = m.positions[i];
		Vec3 B = m.positions[j];
		Vec3 C = m.positions[k];

		/* We must convert float to double as mass() expect Vec3d*/
		Vec3d AB = {(double)B[0] - (double)A[0], (double)B[1] - (double)A[1], (double)B[2] - (double)A[2]};
		Vec3d AC = {(double)B[0] - (double)C[0], (double)B[1] - (double)C[1], (double)B[2] - (double)C[2]};

		double tri_contribution[2];
		mass(AB, AC, tri_contribution);
		M(i,i) += tri_contribution[0];
		M(j,j) += tri_contribution[0];
		M(j,j) += tri_contribution[0];
		M(MAX(i,j), MIN(i,j)) += tri_contribution[1];
		M(MAX(j,k), MIN(j,k)) += tri_contribution[1];
		M(MAX(k,i), MIN(k,i)) += tri_contribution[1];
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

	/* Since we can compute the coefficients for each triangle, we will build M by
	by additionning the contribution of each triangle */

	for (size_t tri_index = 0; tri_index<tri_count; tri_index+=3) { // Triangle ABC d'indices i,j,k
		uint32_t i = m.indices[tri_index];
		uint32_t j = m.indices[tri_index+1];
		uint32_t k = m.indices[tri_index+2];

		Vec3 A = m.positions[i];
		Vec3 B = m.positions[j];
		Vec3 C = m.positions[k];

		/* We must convert float to double as mass() expect Vec3d*/
		Vec3d AB = {(double)B[0] - (double)A[0], (double)B[1] - (double)A[1], (double)B[2] - (double)A[2]};
		Vec3d AC = {(double)B[0] - (double)C[0], (double)B[1] - (double)C[1], (double)B[2] - (double)C[2]};

		double tri_contribution[6];
		mass(AB, AC, tri_contribution);
		S(i,i) += tri_contribution[0];
		S(j,j) += tri_contribution[1];
		S(j,j) += tri_contribution[2];
		S(MAX(i,j), MIN(i,j)) += tri_contribution[3];
		S(MAX(j,k), MIN(j,k)) += tri_contribution[4];
		S(MAX(k,i), MIN(k,i)) += tri_contribution[5];
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
