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

/* CSRMatrix variants */

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
	/* Your implementation goes here.
	 * Use a VTAdjacency structure (see include/matrix/adjacency.h)
	 */
	size_t nv = m.vertex_count();
	size_t nt = m.triangle_count();

	P.symmetric = true;
	P.rows = P.cols = nv;
	P.row_start.resize(nv + 1);

	VTAdjacency vt(m);

	// Majorant du nombre de coefficients non nuls
	const size_t capacity = 3 * nt + nv;
	P.col.resize(capacity);

	size_t write_pos = 0;

	for (size_t i = 0; i < nv; ++i) {
		P.row_start[i] = write_pos;

		uint32_t *row_cols = &P.col[write_pos];
		size_t row_nnz = 0;

		const size_t begin = vt.offset[i];
		const size_t end = begin + vt.degree[i];

		for (size_t k = begin; k < end; ++k) {
			const uint32_t j1 = vt.vtri[k].next;
			const uint32_t j2 = vt.vtri[k].prev;

			if (j1 < i && !find(j1, row_cols, row_nnz)) {
				P.col[write_pos++] = j1;
				++row_nnz;
			}

			if (j2 < i && !find(j2, row_cols, row_nnz)) {
				P.col[write_pos++] = j2;
				++row_nnz;
			}
		}

		// Diagonale
		P.col[write_pos++] = static_cast<uint32_t>(i);
	}

	P.row_start[nv] = write_pos;
	P.col.resize(write_pos);
	P.col.shrink_to_fit();

	// Tri croissant dans chaque ligne
	for (size_t i = 0; i < nv; ++i) {
		uint32_t *row = &P.col[P.row_start[i]];
		const size_t count = P.row_start[i + 1] - P.row_start[i];

		for (size_t k = 1; k < count; ++k) {
			uint32_t val = row[k];
			size_t j = k;

			while (j > 0 && row[j - 1] > val) {
				row[j] = row[j - 1];
				--j;
			}
			row[j] = val;
		}
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
	const TArray<uint32_t> &indices = m.indices;

	for (size_t t = 0; t < tri_count; ++t) {

		const uint32_t ia = indices[3 * t];
		const uint32_t ib = indices[3 * t + 1];
		const uint32_t ic = indices[3 * t + 2];

		const Vec3f &Pa = m.positions[ia];
		const Vec3f &Pb = m.positions[ib];
		const Vec3f &Pc = m.positions[ic];

		Vec3d AB{
			static_cast<double>(Pb[0]) - static_cast<double>(Pa[0]),
			static_cast<double>(Pb[1]) - static_cast<double>(Pa[1]),
			static_cast<double>(Pb[2]) - static_cast<double>(Pa[2])
		};

		Vec3d AC{
			static_cast<double>(Pc[0]) - static_cast<double>(Pa[0]),
			static_cast<double>(Pc[1]) - static_cast<double>(Pa[1]),
			static_cast<double>(Pc[2]) - static_cast<double>(Pa[2])
		};

		double local_M[2];
		mass(AB, AC, local_M);

		M(ia, ia) += local_M[0];
		M(ib, ib) += local_M[0];
		M(ic, ic) += local_M[0];

		const uint32_t ab_max = ia > ib ? ia : ib;
		const uint32_t ab_min = ia > ib ? ib : ia;

		const uint32_t bc_max = ib > ic ? ib : ic;
		const uint32_t bc_min = ib > ic ? ic : ib;

		const uint32_t ca_max = ic > ia ? ic : ia;
		const uint32_t ca_min = ic > ia ? ia : ic;

		M(ab_max, ab_min) += local_M[1];
		M(bc_max, bc_min) += local_M[1];
		M(ca_max, ca_min) += local_M[1];
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

	const TArray<uint32_t> &indices = m.indices;

	for (size_t t = 0; t < tri_count; ++t) {

		const uint32_t ia = indices[3 * t];
		const uint32_t ib = indices[3 * t + 1];
		const uint32_t ic = indices[3 * t + 2];

		const Vec3f &Pa = m.positions[ia];
		const Vec3f &Pb = m.positions[ib];
		const Vec3f &Pc = m.positions[ic];

		Vec3d AB{
			static_cast<double>(Pb[0]) - static_cast<double>(Pa[0]),
			static_cast<double>(Pb[1]) - static_cast<double>(Pa[1]),
			static_cast<double>(Pb[2]) - static_cast<double>(Pa[2])
		};

		Vec3d AC{
			static_cast<double>(Pc[0]) - static_cast<double>(Pa[0]),
			static_cast<double>(Pc[1]) - static_cast<double>(Pa[1]),
			static_cast<double>(Pc[2]) - static_cast<double>(Pa[2])
		};

		double local_S[6];
		stiffness(AB, AC, local_S);

		S(ia, ia) += local_S[0];
		S(ib, ib) += local_S[1];
		S(ic, ic) += local_S[2];

		const uint32_t ab_max = ia > ib ? ia : ib;
		const uint32_t ab_min = ia > ib ? ib : ia;

		const uint32_t bc_max = ib > ic ? ib : ic;
		const uint32_t bc_min = ib > ic ? ic : ib;

		const uint32_t ca_max = ic > ia ? ic : ia;
		const uint32_t ca_min = ic > ia ? ia : ic;

		S(ab_max, ab_min) += local_S[3];
		S(bc_max, bc_min) += local_S[4];
		S(ca_max, ca_min) += local_S[5];
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
