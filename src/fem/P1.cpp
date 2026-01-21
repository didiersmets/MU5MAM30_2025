#include <stdio.h>
#include <string.h>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass_P1.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P1.h"
#include "hash_table.h"
#include "hash.h"

#include <unordered_set>

using namespace std;

/* Neighbors structure for searching ends, behave as std::pair */
struct Neigh
{
	uint32_t first;
	uint32_t second;
	bool operator==(const Neigh &other) const
	{
		return first == other.first && second == other.second;
	}
};

/* Neighbors structure hasher */
struct NeighHasher
{
	static constexpr uint32_t empty_int = ~uint32_t(0);
	static constexpr Neigh empty_key{empty_int, empty_int};

	size_t hash(Neigh nei) const
	{
		uint32_t hash = 0;
		hash = murmur2_32(hash, nei.first);
		hash = murmur2_32(hash, nei.second);
		return hash;
	}

	bool is_empty(Neigh nei) const
	{
		return nei.first == empty_int && nei.second == empty_int;
	}

	bool is_equal(Neigh nei_1, Neigh nei_2) const
	{
		return nei_1.first == nei_2.first && nei_1.second == nei_2.second;
	}
};

/* CSRMatrix variants */
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	/* Adjacency mapping */
	VTAdjacency vt_adj(m);

	size_t nv = m.vertex_count();

	P.symmetric = true;
	P.rows = nv;
	P.cols = nv;

	size_t nnz_non_symmetric = nv;
	for (size_t k = 0; k < nv; k++)
		nnz_non_symmetric += vt_adj.degree[k];

	/* Build row_start and col */
	P.row_start.resize(nv + 1);
	P.col.resize(nnz_non_symmetric);

	NeighHasher hasher{};
	HashTable<Neigh, uint32_t, NeighHasher> seen(nnz_non_symmetric, hasher);

	uint32_t *current_key;
	Neigh current_neigh;
	uint32_t dummy = 0;
	size_t nnz = 0;

	for (size_t k = 0; k < nv; k++)
	{
		P.row_start[k] = nnz;
		for (size_t j = vt_adj.offset[k]; j < vt_adj.offset[k] + vt_adj.degree[k]; j++)
		{
			/* The first neighbor is the current vertex */
			current_neigh.first = k;

			/* The second neighbor is either the vertex, the next or the previous */
			current_neigh.second = k;
			current_key = seen.get(current_neigh);

			if (!current_key)
			{
				P.col[nnz++] = current_neigh.second;
				seen.set_at(current_neigh, dummy);
			}

			current_neigh.second = vt_adj.vtri[j].next;
			current_key = seen.get(current_neigh);

			if (current_neigh.second <= k)
				if (!current_key)
				{
					P.col[nnz++] = current_neigh.second;
					seen.set_at(current_neigh, dummy);
				}

			current_neigh.second = vt_adj.vtri[j].prev;
			current_key = seen.get(current_neigh);

			if (current_neigh.second <= k)
				if (!current_key)
				{
					P.col[nnz++] = current_neigh.second;
					seen.set_at(current_neigh, dummy);
				}
		}
	}
	P.row_start[nv] = nnz;
	P.nnz = nnz;
	P.col.resize(nnz);

	/* Reorder col so that it is consistent with the matrix point of view */
	uint32_t temp;
	for (size_t i = 0; i < nv; i++)
	{
		for (size_t j = 1; j < P.row_start[i + 1] - P.row_start[i]; j++)
			for (size_t k = P.row_start[i]; k < P.row_start[i + 1] - j; k++)
			{
				if (P.col[k + 1] < P.col[k])
				{
					temp = P.col[k];
					P.col[k] = P.col[k + 1];
					P.col[k + 1] = temp;
				}
			}
	}
}

void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
	size_t vtx_count = m.vertex_count();
	assert(P.row_start.size == vtx_count + 1);

	M.symmetric = P.symmetric;
	M.rows = M.cols = vtx_count;
	M.nnz = P.col.size;
	M.row_start = P.row_start.data;
	M.col = P.col.data;
	M.data.resize(M.nnz);
	for (size_t i = 0; i < M.nnz; ++i)
	{
		M.data[i] = 0.0;
	}

	double *mass_matrix = (double *)malloc(3 * 3 * sizeof(double));

	/* Build the global mass matrix */
	size_t nt = m.triangle_count();
	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = m.indices[3 * tri];
		uint32_t b = m.indices[3 * tri + 1];
		uint32_t c = m.indices[3 * tri + 2];

		Vec3 A = m.positions[a];
		Vec3 B = m.positions[b];
		Vec3 C = m.positions[c];

		Vec3d AB = {(double)B[0] - (double)A[0],
					(double)B[1] - (double)A[1],
					(double)B[2] - (double)A[2]};
		Vec3d AC = {(double)C[0] - (double)A[0],
					(double)C[1] - (double)A[1],
					(double)C[2] - (double)A[2]};

		mass_P1(AB, AC, mass_matrix);

		M(a, a) += mass_matrix[0];
		M(b, b) += mass_matrix[4];
		M(c, c) += mass_matrix[8];

		if (b < a)
			M(a, b) += mass_matrix[1];
		else
			M(b, a) += mass_matrix[1];

		if (c < a)
			M(a, c) += mass_matrix[2];
		else
			M(c, a) += mass_matrix[2];

		if (c < b)
			M(b, c) += mass_matrix[5];
		else
			M(c, b) += mass_matrix[5];
	}
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	size_t vtx_count = m.vertex_count();
	assert(P.row_start.size == vtx_count + 1);

	S.symmetric = P.symmetric;
	S.rows = S.cols = vtx_count;
	S.nnz = P.col.size;
	S.row_start = P.row_start.data;
	S.col = P.col.data;
	S.data.resize(S.nnz);
	for (size_t i = 0; i < S.nnz; ++i)
	{
		S.data[i] = 0.0;
	}

	double *stiffness_matrix = (double *)malloc(6 * sizeof(double));

	/* Build the global stiffness matrix */
	size_t nt = m.triangle_count();
	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = m.indices[3 * tri];
		uint32_t b = m.indices[3 * tri + 1];
		uint32_t c = m.indices[3 * tri + 2];

		Vec3 A = m.positions[a];
		Vec3 B = m.positions[b];
		Vec3 C = m.positions[c];

		Vec3d AB = {(double)B[0] - (double)A[0],
					(double)B[1] - (double)A[1],
					(double)B[2] - (double)A[2]};
		Vec3d AC = {(double)C[0] - (double)A[0],
					(double)C[1] - (double)A[1],
					(double)C[2] - (double)A[2]};

		stiffness_P1(AB, AC, stiffness_matrix);

		S(a, a) += stiffness_matrix[0];
		S(b, b) += stiffness_matrix[1];
		S(c, c) += stiffness_matrix[2];

		if (b < a)
			S(a, b) += stiffness_matrix[3];
		else
			S(b, a) += stiffness_matrix[3];

		if (c < a)
			S(a, c) += stiffness_matrix[5];
		else
			S(c, a) += stiffness_matrix[5];

		if (c < b)
			S(b, c) += stiffness_matrix[4];
		else
			S(c, b) += stiffness_matrix[4];
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
	for (size_t t = 0; t < tri_count; ++t)
	{
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = {(double)B[0] - (double)A[0],
					(double)B[1] - (double)A[1],
					(double)B[2] - (double)A[2]};
		Vec3d AC = {(double)C[0] - (double)A[0],
					(double)C[1] - (double)A[1],
					(double)C[2] - (double)A[2]};
		double Mloc[2];
		mass_P1(AB, AC, Mloc);
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
	for (size_t t = 0; t < tri_count; ++t)
	{
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = {(double)B[0] - (double)A[0],
					(double)B[1] - (double)A[1],
					(double)B[2] - (double)A[2]};
		Vec3d AC = {(double)C[0] - (double)A[0],
					(double)C[1] - (double)A[1],
					(double)C[2] - (double)A[2]};
		double Sloc[6];
		stiffness_P1(AB, AC, Sloc);
		S.diag[a] += Sloc[0];
		S.diag[b] += Sloc[1];
		S.diag[c] += Sloc[2];
		S.off_diag[3 * t + 0] = Sloc[3];
		S.off_diag[3 * t + 1] = Sloc[4];
		S.off_diag[3 * t + 2] = Sloc[5];
	}
}
