#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P2.h"
#include "edge.h"
#include "adjacency.h"
#include "mass.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness.h"
#include "math_utils.h"


// Auxiliary function
static bool is_in_first_idx(uint32_t target, uint32_t *array, size_t nb_idx_check)
{
	for (size_t i=0; i<nb_idx_check; i++) {
		if (array[i]==target) {
			return true;
		}
	}
	return false;
}

// Fonction lambda pour retrouver l'indice P2 (DDL) d'un point milieu par recherche dichotomique
auto get_midpoint_id(uint32_t a, uint32_t b, std::vector<Edge> &out_edges, size_t offset) {
    Edge e = {std::min(a, b), std::max(a, b)};
    auto it = std::lower_bound(out_edges.begin(), out_edges.end(), e);
    return offset + std::distance(out_edges.begin(), it);
};

void build_P2_CSRPattern(const Mesh &m, CSRPattern &P, std::vector<Edge> &out_edges)
{
    out_edges.reserve(m.triangle_count() * 3);

    for (size_t t = 0; t < m.index_count(); t+=3) {

        uint32_t i = m.indices[t];
        uint32_t j = m.indices[t+1];
        uint32_t k = m.indices[t+2];

        out_edges.push_back({MIN(i, j), MAX(i, j)});
        out_edges.push_back({MIN(j, k), MAX(j, k)});
        out_edges.push_back({MIN(k, i), MAX(k, i)});
    }

    // Tri et suppression des doublons (O(N log N))
    std::sort(out_edges.begin(), out_edges.end());
    out_edges.erase(std::unique(out_edges.begin(), out_edges.end()), out_edges.end());

    VTAdjacency vtadj(m);
	P.symmetric = true;
	P.rows = m.vertex_count() + out_edges.size();
	P.cols = m.vertex_count() + out_edges.size();


	/* Note that each vertex index i correspond to a line in the sparse matrix 
	   therefore we build the pattern line by line */
	P.row_start.resize(P.rows+1);
	P.col.resize(15*m.triangle_count() + P.rows); // upper estimate

	size_t nnz = 0;
    // Iterate on vertex first
	for (uint32_t i = 0; i<m.vertex_count(); i++) {
		P.row_start[i] = nnz; // The i line begins at this nnz-th value
		size_t line_nnz = 0;
		uint32_t *line_start = &P.col[nnz];

		/* iterate on every triangle connected to vertex i */
		uint32_t tri_start = vtadj.offset[i]; // first triangle in VTAdjency structure
		uint32_t tri_stop = tri_start + vtadj.degree[i]; // next to last triangle

		for (size_t tri_index = tri_start; tri_index<tri_stop; tri_index++) { // triangle ijk
			uint32_t j = vtadj.vtri[tri_index].next;
			uint32_t k = vtadj.vtri[tri_index].prev;
            
            // mid-points of ijk triangle
            uint32_t ij = get_midpoint_id(i, j, out_edges, m.vertex_count());
            uint32_t ik = get_midpoint_id(i, k, out_edges, m.vertex_count());
            uint32_t jk = get_midpoint_id(j, k, out_edges, m.vertex_count());

            uint32_t dof[5] = {j, k, ij, ik, jk}; // other dof of the triangle

            for (size_t t=0; t<5; t++) {
                if (dof[t]<i && !is_in_first_idx(dof[t], line_start, line_nnz)) { // Check if the dof was not already encountered
                    line_nnz++;
                    P.col[nnz++] = dof[t];
			    }
            }
		}
		P.col[nnz++] = i; // i is always connected to himself (diag term)
	}
    // Iterate on the mid-points

    for (size_t t=0; t<out_edges.size(); t++) {
        uint32_t mid_point_id = m.vertex_count() + t;
        P.row_start[mid_point_id] = nnz;
        size_t line_nnz = 0;
        uint32_t *line_start = &P.col[nnz];

        // Get the two endpoint generating the midpoint
        uint32_t i = out_edges[t].first;
        uint32_t j = out_edges[t].second;

        // We only check adjacency on i, if j appears we get the third vertex
        uint32_t tri_start = vtadj.offset[i];
		uint32_t tri_stop = tri_start + vtadj.degree[i];

        for (size_t tri_index=tri_start; tri_index<tri_stop; tri_index++) {
            uint32_t next_vertex = vtadj.vtri[tri_index].next;
            uint32_t prev_vertex = vtadj.vtri[tri_index].prev;

            if (next_vertex != j && prev_vertex !=j ) {continue;} // check if j is in the triangle
            uint32_t k = (next_vertex == j) ? prev_vertex : next_vertex;

            // get the other mid-points
            uint32_t ik = get_midpoint_id(i, k, out_edges, m.vertex_count());
            uint32_t jk = get_midpoint_id(j, k, out_edges, m.vertex_count());

            uint32_t dof[5] = {i, j, k, ik, jk};

            for (size_t idx=0; idx<5; idx++) {
                if (dof[idx]<mid_point_id && !is_in_first_idx(dof[idx], line_start, line_nnz)) {
                    line_nnz++;
                    P.col[nnz++] = dof[idx];
                }
            }
        }
        P.col[nnz++] = mid_point_id;
    }

	P.row_start[P.rows] = nnz;
	P.nnz = nnz;
	P.col.resize(nnz); // last size was upper estimate

	/* each col needs to be sorted in the final pattern*/
	for (size_t i = 0; i < P.rows; i++) {
		size_t line_start = P.row_start[i];
		size_t line_end = P.row_start[i+1];

		std::sort(P.col.data + line_start, P.col.data + line_end); // use a O(N*log(N)) sorting algorithm
	}
}


void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M, std::vector<Edge> &out_edges)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	assert(P.row_start.size == vtx_count + out_edges.size() + 1);

	M.symmetric = true;
	M.rows = M.cols = vtx_count + out_edges.size();
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
	for (size_t tri_index = 0; tri_index<tri_count; tri_index++) { // Triangle ABC d'indices i,j,k
		uint32_t i = m.indices[3*tri_index];
		uint32_t j = m.indices[3*tri_index+1];
		uint32_t k = m.indices[3*tri_index+2];

        uint32_t ij = get_midpoint_id(i, j, out_edges, m.vertex_count());
        uint32_t jk = get_midpoint_id(j, k, out_edges, m.vertex_count());
        uint32_t ik = get_midpoint_id(i, k, out_edges, m.vertex_count());

		Vec3 A = m.positions[i];
		Vec3 B = m.positions[j];
		Vec3 C = m.positions[k];

		/* We must convert float to double as mass() expect Vec3d */
		Vec3d AB = {(double)B[0] - (double)A[0], (double)B[1] - (double)A[1], (double)B[2] - (double)A[2]};
		Vec3d AC = {(double)C[0] - (double)A[0], (double)C[1] - (double)A[1], (double)C[2] - (double)A[2]};

		double tri_contribution[5];
		mass_P2(AB, AC, tri_contribution);
		M(i,i) += tri_contribution[0];
		M(j,j) += tri_contribution[0];
		M(k,k) += tri_contribution[0];
		M(MAX(i,j), MIN(i,j)) += tri_contribution[1];
		M(MAX(j,k), MIN(j,k)) += tri_contribution[1];
		M(MAX(k,i), MIN(k,i)) += tri_contribution[1];
        M(ij,ij) += tri_contribution[2];
        M(jk,jk) += tri_contribution[2];
        M(ik,ik) += tri_contribution[2];
        M(MAX(ij,jk), MIN(ij,jk)) += tri_contribution[3];
        M(MAX(jk,ik), MIN(jk,ik)) += tri_contribution[3];
        M(MAX(ij,ik), MIN(ij,ik)) += tri_contribution[3];
        M(MAX(i,jk), MIN(i,jk)) += tri_contribution[4];
        M(MAX(j,ik), MIN(j,ik)) += tri_contribution[4];
        M(MAX(k,ij), MIN(k,ij)) += tri_contribution[4];
	}
}

void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S, std::vector<Edge> &out_edges)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	assert(P.row_start.size == vtx_count + out_edges.size() + 1);

	S.symmetric = true;
	S.rows = S.cols = vtx_count + out_edges.size();
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

	for (size_t tri_index = 0; tri_index<tri_count; tri_index++) { // Triangle ABC d'indices i,j,k
		uint32_t i = m.indices[3*tri_index];
		uint32_t j = m.indices[3*tri_index+1];
		uint32_t k = m.indices[3*tri_index+2];

        uint32_t ij = get_midpoint_id(i, j, out_edges, m.vertex_count());
        uint32_t jk = get_midpoint_id(j, k, out_edges, m.vertex_count());
        uint32_t ik = get_midpoint_id(i, k, out_edges, m.vertex_count());

		Vec3 A = m.positions[i];
		Vec3 B = m.positions[j];
		Vec3 C = m.positions[k];

		/* We must convert float to double as mass() expect Vec3d*/
		Vec3d AB = {(double)B[0] - (double)A[0], (double)B[1] - (double)A[1], (double)B[2] - (double)A[2]};
		Vec3d AC = {(double)C[0] - (double)A[0], (double)C[1] - (double)A[1], (double)C[2] - (double)A[2]};

		double tri_contribution[15];
		stiffness_P2(AB, AC, tri_contribution);
		S(i,i) += tri_contribution[0];
		S(j,j) += tri_contribution[1];
		S(k,k) += tri_contribution[2];
        S(ij,ij) += tri_contribution[3];
        S(jk,jk) += tri_contribution[4];
        S(ik,ik) += tri_contribution[5];
		S(MAX(i,j), MIN(i,j)) += tri_contribution[6];
		S(MAX(j,k), MIN(j,k)) += tri_contribution[7];
		S(MAX(k,i), MIN(k,i)) += tri_contribution[8];
        S(MAX(ij,jk), MIN(ij,jk)) += tri_contribution[9];
        S(MAX(jk,ik), MIN(jk,ik)) += tri_contribution[10];
        S(MAX(ij,ik), MIN(ij,ik)) += tri_contribution[11];
        S(MAX(i,ij), MIN(i,ij)) += tri_contribution[12];
        S(MAX(j,ij), MIN(j,ij)) += tri_contribution[12];
        S(MAX(j,jk), MIN(j,jk)) += tri_contribution[13];
        S(MAX(k,jk), MIN(k,jk)) += tri_contribution[13];
        S(MAX(i,ik), MIN(i,ik)) += tri_contribution[14];
        S(MAX(k,ik), MIN(k,ik)) += tri_contribution[14];
	}
}