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

#define INCLUDE_DIAGONAL 1

/* CSRMatrix variants */

// void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
// {
// 	VTAdjacency vta(m);
// 	P.nnz = m.triangle_count() * 3 + INCLUDE_DIAGONAL * m.vertex_count();
// 	P.row_start.resize(m.vertex_count() + 1);
// 	P.col.resize(P.nnz);
// 	P.rows = m.vertex_count();
// 	P.cols = m.vertex_count();

// 	// Setting up row_start
// 	P.row_start[0] = 0;
// 	for ( size_t i = 0; i < P.rows; i++ ) {
// 		P.row_start[i] = vta.offset[i] + INCLUDE_DIAGONAL * i;  // + 1 if I count also the diagonal
// 	}
// 	P.row_start[m.vertex_count()] = P.row_start[m.vertex_count()-1] + vta.degree[m.vertex_count()-1] + INCLUDE_DIAGONAL;
	

// 	// Setting up cols
// 	for ( size_t i = 0; i < P.rows; i++ ) {
// 		// Scrolling through each neighbour
// 		for ( size_t j = 0; j < vta.degree[i]; j++ ) {
// 			P.col[P.row_start[i]+j] = vta.vtri[vta.offset[i]+j].next;			// corrected [j] -> [vta.offset[i]]
// 		}
// 		if ( INCLUDE_DIAGONAL ) {
// 			P.col[P.row_start[i]+vta.degree[i]] = i;
// 		}
// 		// Reordering neighbourhood
// 		for (size_t j = 0; j < vta.degree[i] + INCLUDE_DIAGONAL; j++) {
// 			size_t idx = j;

// 			for (size_t k = j + 1; k < vta.degree[i] + INCLUDE_DIAGONAL; k++) {
// 				if (P.col[P.row_start[i] + k] <
// 					P.col[P.row_start[i] + idx]) {
// 					idx = k;   // salva posizione del minimo
// 				}
// 			}

// 			if (idx != j) {
// 				std::swap(
// 					P.col[P.row_start[i] + j],
// 					P.col[P.row_start[i] + idx]
// 				);
// 			}
// 		}
// 	}
// }

static bool find(uint32_t x, uint32_t *start, size_t count)
{
	for (size_t i = 0; i < count; ++i) {
		if (start[i] == x)
			return true;
	}
	return false;
}

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();

	P.symmetric = true;
	P.rows = P.cols = vtx_count;
	P.row_start.resize(vtx_count + 1);

	VTAdjacency adj(m);

	/* Upper bound on the number of edges a->b with a <= b */
	size_t max_nnz = 3 * tri_count + vtx_count;
	P.col.resize(max_nnz);

	/* Fill P.row_start and P.col (not yet ordered) */
	size_t nnz = 0;
	for (size_t a = 0; a < vtx_count; ++a) {
		P.row_start[a] = nnz;
		uint32_t *start = &P.col[nnz];
		size_t nnz_loc = 0;
		uint32_t kstart = adj.offset[a];
		uint32_t kstop = kstart + adj.degree[a];
		for (size_t k = kstart; k < kstop; ++k) {
			uint32_t b = adj.vtri[k].next;
			uint32_t c = adj.vtri[k].prev;
			if (b < a && !find(b, start, nnz_loc)) {
				P.col[nnz++] = b;
				nnz_loc++;
			}
			if (c < a && !find(c, start, nnz_loc)) {
				P.col[nnz++] = c;
				nnz_loc++;
			}
		}
		P.col[nnz++] = a;
	}
	P.row_start[vtx_count] = nnz;
	P.col.resize(nnz);
	P.col.shrink_to_fit();

	/* Reorder each "line" of P.col in increasing column order */
	for (size_t a = 0; a < vtx_count; ++a) {
		uint32_t *__restrict to_sort = &P.col[P.row_start[a]];
		size_t count = P.row_start[a + 1] - P.row_start[a];
		/* Insertion sort (small degree vtx) */
		for (size_t k = 1; k < count; ++k) {
			size_t j = k;
			while (j && (to_sort[j - 1] > to_sort[j])) {
				uint32_t tmp = to_sort[j - 1];
				to_sort[j - 1] = to_sort[j];
				to_sort[j] = tmp;
				j--;
			}
		}
	}
}


void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
	// A "relationship" between two vertices is an edge, which means it
	// belongs to 2 triangles. It means that if I iterate over all the 
	// triangles and update the matrix, I would add twice each time the relationship
	// Then I have to half it.
	// The only time I wouldn't add twice an edge is if it is a border edge,
	// which cannot happen on a sphere

	size_t vtx_count = m.vertex_count();
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


	float* M_loc;

	TVec3<float> AB;
	TVec3<float> AC;
	size_t tria_indices[3] = {0,0,0};  // {A_index, B_index, C_index}

	M_loc = (float *) std::malloc(9 * sizeof(float));

	for ( size_t tria_index = 0; tria_index < m.index_count(); tria_index += 3 ) {
		tria_indices[0] = m.indices[tria_index];    // A_index
		tria_indices[1] = m.indices[tria_index+1];  // B_index
		tria_indices[2] = m.indices[tria_index+2];  // C_index

		AB = m.positions[tria_indices[1]] - m.positions[tria_indices[0]];
		AC = m.positions[tria_indices[2]] - m.positions[tria_indices[0]];

		if ( AB.x == 0 && AB.y == 0 && AB.z == 0 ) {
			printf("AB is zero. A=B=(%f,%f,%f),\n", m.positions[tria_indices[1]].x, m.positions[tria_indices[1]].y, m.positions[tria_indices[1]].z);
		}
		mass<float>(AB, AC, M_loc);

		// M[i_glob][j_glob] += M_loc[i_loc][j_loc]
		for ( int i_loc = 0; i_loc < 3; i_loc++ ) {
			for ( int j_loc = 0; j_loc < 3; j_loc++ ) {
				if ( tria_indices[i_loc] >= tria_indices[j_loc] ) {
					M(tria_indices[i_loc], tria_indices[j_loc]) += M_loc[i_loc*3+j_loc];
				}
			}
		}
	}

	free(M_loc);
	// size_t vtx_count = m.vertex_count();
	// size_t tri_count = m.triangle_count();
	// assert(P.row_start.size == vtx_count + 1);

	// M.symmetric = true;
	// M.rows = M.cols = vtx_count;
	// M.nnz = P.col.size;
	// M.row_start = P.row_start.data;
	// M.col = P.col.data;
	// M.data.resize(M.nnz);
	// for (size_t i = 0; i < M.nnz; ++i) {
	// 	M.data[i] = 0.0;
	// }

	// const TArray<uint32_t> &idx = m.indices;
	// for (size_t t = 0; t < tri_count; ++t) {
	// 	uint32_t a = idx[3 * t + 0];
	// 	uint32_t b = idx[3 * t + 1];
	// 	uint32_t c = idx[3 * t + 2];
	// 	Vec3f A = m.positions[a];
	// 	Vec3f B = m.positions[b];
	// 	Vec3f C = m.positions[c];
	// 	Vec3d AB = { (double)B[0] - (double)A[0],
	// 		     (double)B[1] - (double)A[1],
	// 		     (double)B[2] - (double)A[2] };
	// 	Vec3d AC = { (double)C[0] - (double)A[0],
	// 		     (double)C[1] - (double)A[1],
	// 		     (double)C[2] - (double)A[2] };
	// 	double Mloc[2];
	// 	mass(AB, AC, Mloc);
	// 	M(a, a) += Mloc[0];
	// 	M(b, b) += Mloc[0];
	// 	M(c, c) += Mloc[0];
	// 	M(a > b ? a : b, a > b ? b : a) += Mloc[1];
	// 	M(b > c ? b : c, b > c ? c : b) += Mloc[1];
	// 	M(c > a ? c : a, c > a ? a : c) += Mloc[1];
	// }
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	// size_t vtx_count = m.vertex_count();
	// assert(P.row_start.size == vtx_count + 1);

	// S.symmetric = true;
	// S.rows = S.cols = vtx_count;
	// S.nnz = P.col.size;
	// S.row_start = P.row_start.data;
	// S.col = P.col.data;
	// S.data.resize(S.nnz);
	// for (size_t i = 0; i < S.nnz; ++i) {
	// 	S.data[i] = 0.0;
	// }

	// float* S_loc;

	// TVec3<float> AB;
	// TVec3<float> AC;
	// size_t tria_indices[3] = {0,0,0};  // {A_index, B_index, C_index}

	// S_loc = (float *) std::malloc(9 * sizeof(float));

	// for ( size_t tria_index = 0; tria_index < m.index_count(); tria_index += 3 ) {
	// 	tria_indices[0] = m.indices[tria_index];    // A_index
	// 	tria_indices[1] = m.indices[tria_index+1];  // B_index
	// 	tria_indices[2] = m.indices[tria_index+2];  // C_index

	// 	AB = m.positions[tria_indices[1]] - m.positions[tria_indices[0]];
	// 	AC = m.positions[tria_indices[2]] - m.positions[tria_indices[0]];

	// 	stiffness<float>(AB, AC, S_loc);

	// 	// M[i_glob][j_glob] += M_loc[i_loc][j_loc]
	// 	for ( int i_loc = 0; i_loc < 3; i_loc++ ) {
	// 		for ( int j_loc = 0; j_loc < 3; j_loc++ ) {
	// 			if ( tria_indices[i_loc] >= tria_indices[j_loc] ) {
	// 				S(tria_indices[i_loc], tria_indices[j_loc]) += S_loc[i_loc*3 + j_loc];
	// 			}
	// 		}
	// 	}
	// }

	// free(S_loc);
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
		S(a, a) += Sloc[0];
		S(b, b) += Sloc[1];
		S(c, c) += Sloc[2];
		S(a > b ? a : b, a > b ? b : a) += Sloc[3];
		S(b > c ? b : c, b > c ? c : b) += Sloc[4];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5];
	}
}

/* FEMatrix variants */
// my implementation

void build_P1_mass_matrix(const Mesh &m, FEMatrix &M) {
	// double* M_loc;

	// TVec3<double> AB;
	// TVec3<double> AC;
	// size_t tria_indices[3] = {0,0,0};  // {A_index, B_index, C_index}

	// M_loc = (double *) std::malloc(9 * sizeof(double));

	// for ( int tria_index = 0; tria_index < m.index_count(); tria_index += 3 ) {
	// 	tria_indices[0] = m.indices[tria_index];    // A_index
	// 	tria_indices[1] = m.indices[tria_index+1];  // B_index
	// 	tria_indices[2] = m.indices[tria_index+2];  // C_index

	// 	AB = m.positions[tria_indices[1]] - m.positions[m.indices[tria_indices[0]]];
	// 	AC = m.positions[tria_indices[2]] - m.positions[m.indices[tria_indices[0]]];

	// 	mass(AB, AC, M_loc);

	// 	// M[i_glob][j_glob] += M_loc[i_loc][j_loc]
	// 	for ( int i_loc = 0; i_loc < 3; i_loc++ ) {
	// 		for ( int j_loc = 0; j_loc < 3; j_loc++ ) {
	// 			M[tria_indices[i_loc]][tria_indices[j_loc]] = M_loc[i_loc*3+j_loc];
	// 		}
	// 	}
	// }
}

void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S) {
	// double* S_loc;

	// TVec3<double> AB;
	// TVec3<double> AC;
	// size_t tria_indices[3] = {0,0,0};  // {A_index, B_index, C_index}

	// M_loc = (double *) std::malloc(9 * sizeof(double));

	// for ( int tria_index = 0; tria_index < m.index_count; tria_index += 3 ) {
	// 	tria_indices[0] = m.indices[tria_index];    // A_index
	// 	tria_indices[1] = m.indices[tria_index+1];  // B_index
	// 	tria_indices[2] = m.indices[tria_index+2];  // C_index

	// 	AB = m.positions[tria_indices[1]] - m.positions[m.indices[tria_indices[0]]];
	// 	AC = m.positions[tria_indices[2]] - m.positions[m.indices[tria_indices[0]]];

	// 	stiffness(AB, AC, S_loc);

	// 	// M[i_glob][j_glob] += M_loc[i_loc][j_loc]
	// 	for ( int i_loc = 0; i_loc < 3; i_loc++ ) {
	// 		for ( int j_loc = 0; j_loc < 3; j_loc++ ) {
	// 			S[tria_indices[i_loc]][tria_indices[j_loc]] = S_loc[i_loc*3 + j_loc];
	// 		}
	// 	}
	// }
}

/* FEMatrix variants */
// prof implementations


// void build_P1_mass_matrix(const Mesh &m, FEMatrix &M)
// {
// 	size_t vtx_count = m.vertex_count();
// 	size_t tri_count = m.triangle_count();

// 	M.fem_type = FEMatrix::P1_cst;
// 	M.m = &m;
// 	M.rows = M.cols = vtx_count;

// 	M.diag.resize(vtx_count);
// 	memset(M.diag.data, 0, vtx_count * sizeof(double));

// 	M.off_diag.resize(tri_count);
// 	const TArray<uint32_t> &idx = m.indices;
// 	for (size_t t = 0; t < tri_count; ++t) {
// 		uint32_t a = idx[3 * t + 0];
// 		uint32_t b = idx[3 * t + 1];
// 		uint32_t c = idx[3 * t + 2];
// 		Vec3f A = m.positions[a];
// 		Vec3f B = m.positions[b];
// 		Vec3f C = m.positions[c];
// 		Vec3d AB = { (double)B[0] - (double)A[0],
// 			     (double)B[1] - (double)A[1],
// 			     (double)B[2] - (double)A[2] };
// 		Vec3d AC = { (double)C[0] - (double)A[0],
// 			     (double)C[1] - (double)A[1],
// 			     (double)C[2] - (double)A[2] };
// 		double Mloc[2];
// 		mass(AB, AC, Mloc);
// 		M.diag[a] += Mloc[0];
// 		M.diag[b] += Mloc[0];
// 		M.diag[c] += Mloc[0];
// 		M.off_diag[t] = Mloc[1];
// 	}
// }

// void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S)
// {
// 	size_t vtx_count = m.vertex_count();
// 	size_t tri_count = m.triangle_count();

// 	S.fem_type = FEMatrix::P1_sym;
// 	S.m = &m;
// 	S.rows = S.cols = vtx_count;

// 	S.diag.resize(vtx_count);
// 	memset(S.diag.data, 0, vtx_count * sizeof(double));

// 	S.off_diag.resize(3 * tri_count);
// 	const TArray<uint32_t> &idx = m.indices;
// 	for (size_t t = 0; t < tri_count; ++t) {
// 		uint32_t a = idx[3 * t + 0];
// 		uint32_t b = idx[3 * t + 1];
// 		uint32_t c = idx[3 * t + 2];
// 		Vec3f A = m.positions[a];
// 		Vec3f B = m.positions[b];
// 		Vec3f C = m.positions[c];
// 		Vec3d AB = { (double)B[0] - (double)A[0],
// 			     (double)B[1] - (double)A[1],
// 			     (double)B[2] - (double)A[2] };
// 		Vec3d AC = { (double)C[0] - (double)A[0],
// 			     (double)C[1] - (double)A[1],
// 			     (double)C[2] - (double)A[2] };
// 		double Sloc[6];
// 		stiffness(AB, AC, Sloc);
// 		S.diag[a] += Sloc[0];
// 		S.diag[b] += Sloc[1];
// 		S.diag[c] += Sloc[2];
// 		S.off_diag[3 * t + 0] = Sloc[3];
// 		S.off_diag[3 * t + 1] = Sloc[4];
// 		S.off_diag[3 * t + 2] = Sloc[5];
// 	}
// }
