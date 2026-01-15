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

#include <unordered_set>
#include <algorithm>

/* CSRMatrix variants */

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
  /* (done) Your implementation goes here.
   * Use a VTAdjacency structure (see include/matrix/adjacency.h)
   */

  size_t n_vtx = m.vertex_count();
  VTAdjacency A(m);

  P.symmetric = true;
  P.rows = n_vtx;
  P.cols = n_vtx;
  P.nnz = A.vtri.size / 2 + n_vtx;
  P.row_start.resize(P.rows + 1);
  P.col.resize(P.nnz);

  uint32_t col_pos = 0;
  P.row_start[0] = 0;
  for (size_t i=0; i<n_vtx; ++i) {

    std::unordered_set<uint32_t> big_neighbs;
    for (uint32_t k=A.offset[i]; k<A.offset[i]+A.degree[i]; ++k) {
      uint32_t j;

      j = A.vtri[k].prev;
      if (i<j)
	big_neighbs.insert(j);

      j = A.vtri[k].next;
      if (i<j)
	big_neighbs.insert(j);
    }
    for (const uint32_t &j : big_neighbs)
	P.col[col_pos++] = j;

    /* Put the i index at then end (needed for CSRMatrix::sum) */
    P.col[col_pos++] = i;

    P.row_start[i+1] = col_pos;
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

  /* (done) Your implementation goes here */

  for (size_t tri=0; tri<tri_count; ++tri) {
    uint32_t v_ids[3] = {m.indices[3*tri],
			 m.indices[3*tri+1],
			 m.indices[3*tri+2]};
    std::sort(std::begin(v_ids), std::end(v_ids));
    uint32_t va = v_ids[0];
    uint32_t vb = v_ids[1];
    uint32_t vc = v_ids[2];

    Vec3f A = m.positions[va];
    Vec3f B = m.positions[vb];
    Vec3f C = m.positions[vc];

    Vec3d AB = { (double)(B[0] - A[0]),
		 (double)(B[1] - A[1]),
		 (double)(B[2] - A[2])};
    Vec3d AC = { (double)(C[0] - A[0]),
		 (double)(C[1] - A[1]),
		 (double)(C[2] - A[2])};

    double M_tri[2];
    mass(AB, AC, M_tri);

    M(va, va) += M_tri[0];
    M(vb, vb) += M_tri[0];
    M(vc, vc) += M_tri[0];
    M(va, vb) += M_tri[1];
    M(va, vc) += M_tri[1];
    M(vb, vc) += M_tri[1];
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

  /* (done) Your implementation goes here */

  for (size_t tri=0; tri<tri_count; ++tri) {
    uint32_t v_ids[3] = {m.indices[3*tri],
			 m.indices[3*tri+1],
			 m.indices[3*tri+2]};
    std::sort(std::begin(v_ids), std::end(v_ids));
    uint32_t va = v_ids[0];
    uint32_t vb = v_ids[1];
    uint32_t vc = v_ids[2];

    Vec3f A = m.positions[va];
    Vec3f B = m.positions[vb];
    Vec3f C = m.positions[vc];

    Vec3d AB = { (double)(B[0] - A[0]),
		 (double)(B[1] - A[1]),
		 (double)(B[2] - A[2])};
    Vec3d AC = { (double)(C[0] - A[0]),
		 (double)(C[1] - A[1]),
		 (double)(C[2] - A[2])};

    double S_tri[6];
    stiffness(AB, AC, S_tri);

    S(va, va) += S_tri[0];
    S(vb, vb) += S_tri[1];
    S(vc, vc) += S_tri[2];
    S(va, vb) += S_tri[3];
    S(vb, vc) += S_tri[4];
    S(va, vc) += S_tri[5];
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
