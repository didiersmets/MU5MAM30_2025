#include <stdio.h>
#include <string.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "adjacency.h"
#include "mass_P2.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P2.h"

#include <set>
#include <algorithm>
#include <map>
#include <unordered_map>

/* CSRMatrix variants */

void build_P2_CSRPattern(Mesh &m,
			 CSRPattern &P)
{
  /*
    Global dofs numbering strategy:
    k in {0,...,n_vtx-1} : corresponds to vertex k
    k in {n_vtx,...,n_vtx+n_edges-1}: apply the following indexing algorithm
    - set k = n_vtx
    - loop over all vertices v from 0 to n_vtx in increasing index order
    - loop over the neighbors u of v in increasing index order
    - if v < u: label the middle of edge (u,v) with index k and then do k++
    This procedure gives exactly the lexicographical order on the edges with edge representation
    (u,v) with u<v.

    Symetric matrix storage : dofs (k,l) such that k >= l
  */

#if USE_P2
  auto &edge2dof = m.edge2dof;
#else
  std::unordered_map<uint64_t, uint32_t> edge2dof;
#endif

  size_t n_vtx = m.vertex_count();
  size_t n_tri = m.triangle_count();
  size_t n_edges = m.index_count() / 2;
  size_t n_dofs = n_vtx + n_edges;
  VTAdjacency A(m);

  P.symmetric = true;
  P.rows = n_dofs;
  P.cols = n_dofs;
  P.row_start.resize(P.rows + 1);
  /*
    P.nnz =
    n_dofs       +  // (i,i) coefficients
    n_edges      +  // (i,j) coeffs for i,j vertices
    2 * n_edges  +  // (i,j) coeffs for i vertex and j edge middle
    4 * n_edges;    // (i,j) coeffs for i,j edges middles
    P.col.resize(P.nnz);
  */

  uint32_t edge_dof = n_vtx;
  P.row_start[0] = 0;
  for (size_t i=0; i<n_vtx; ++i) {

    /* Get all neighbors of i and order them in increasing order */
    std::set<uint32_t> neighbs;
    for (uint32_t k=A.offset[i]; k<A.offset[i]+A.degree[i]; ++k) {
      neighbs.insert(A.vtri[k].prev);
      neighbs.insert(A.vtri[k].next);
    }

    /* Number the edges */
    for (const uint32_t &j : neighbs) {
      if (i<j)
	edge2dof.insert({pack(i, j), edge_dof++});
    }

    /*
      Add the coefficients (we know that edge middle dofs have bigger indexes than
      vertex dofs so dont add edge middle <-> vertex interactions here)
    */
    for (const uint32_t &j : neighbs) {
      if (j >= i)
	break;
      P.col.push_back(j);
    }

    /* Put the i index at the end (needed for CSRMatrix::sum) */
    P.col.push_back(i);

    P.row_start[i+1] = P.col.size;
  }
  assert(edge_dof == n_dofs);

  /*
    At this point we have added all vertex <-> vertex coeffcients and we have numbered all
    edge middle points. We need to add edge middle <-> vertex and ege middle <-> edge middle
    coefficients.
  */

  /*
    Now we loop over all triangles and add the edges middle points to other dofs iteractions
    Since we don't know a priori how many non zero coeffs will be stored on each line for
    each edge middle we first put the interactions in a multimap
  */
  std::multimap<uint64_t, uint32_t> e2lowerdof;
  for (size_t tri=0; tri<n_tri; ++tri) {
    uint32_t v_ids[3] = {m.indices[3*tri],
			 m.indices[3*tri+1],
			 m.indices[3*tri+2]};
    std::sort(std::begin(v_ids), std::end(v_ids));
    uint32_t va = v_ids[0];
    uint32_t vb = v_ids[1];
    uint32_t vc = v_ids[2];

    uint32_t eab = edge2dof.find(pack(va, vb))->second;
    uint32_t eac = edge2dof.find(pack(va, vc))->second;
    uint32_t ebc = edge2dof.find(pack(vb, vc))->second;

    /* This is ordered in increasing dof id */
    uint32_t dof_ids[6] = {va, vb, vc, eab, eac, ebc};

    for (uint32_t e_loc_id=3; e_loc_id<6; ++e_loc_id) {
      uint32_t e_dof = dof_ids[e_loc_id];
      for (uint32_t k=0; k<e_loc_id; ++k) {
	assert(e_dof > dof_ids[k]);
	e2lowerdof.insert({e_dof, dof_ids[k]});
      }
    }
  }

  /* Now we fill up P with the corresponding edge middle <-> everything else interactions */
  uint64_t curr_e = n_vtx;
  for (const auto &p : e2lowerdof) {
    if (p.first != curr_e) {
      assert(p.first == curr_e+1);
      /* Put the edge middle to itself interaction at the end (needed for CSRMatrix::sum) */
      P.col.push_back(curr_e);
      P.row_start[++curr_e] = P.col.size;
      assert(curr_e == p.first);
    }
    P.col.push_back(p.second);
  }
  P.col.push_back(curr_e);
  P.row_start[++curr_e] = P.col.size;

  P.nnz = P.col.size;
  assert(curr_e == n_dofs);
}

void build_P2_mass_matrix(const Mesh &m,
			  const CSRPattern &P,
			  CSRMatrix &M)
{
#if USE_P2
  auto &edge2dof = m.edge2dof;
#else
  std::unordered_map<uint64_t, uint32_t> edge2dof;
#endif

  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();
  size_t n_edges = 3 * tri_count / 2;
  assert(P.row_start.size == vtx_count + n_edges + 1);

  M.symmetric = true;
  M.rows = P.rows;
  M.cols = P.cols;
  M.nnz = P.nnz;;
  M.row_start = P.row_start.data;
  M.col = P.col.data;
  M.data.resize(M.nnz);
  for (size_t i = 0; i < M.nnz; ++i) {
    M.data[i] = 0.0;
  }

  for (size_t tri=0; tri<tri_count; ++tri) {
    uint32_t v_ids[3] = {m.indices[3*tri],
			 m.indices[3*tri+1],
			 m.indices[3*tri+2]};
    std::sort(std::begin(v_ids), std::end(v_ids));
    uint32_t va = v_ids[0];
    uint32_t vb = v_ids[1];
    uint32_t vc = v_ids[2];

    uint32_t eab = edge2dof.find(pack(va, vb))->second;
    uint32_t eac = edge2dof.find(pack(va, vc))->second;
    uint32_t ebc = edge2dof.find(pack(vb, vc))->second;

    uint32_t dof_ids[6] = {va, vb, vc, eab, eac, ebc};

    Vec3f A = m.positions[va];
    Vec3f B = m.positions[vb];
    Vec3f C = m.positions[vc];

    Vec3d AB = { (double)(B[0] - A[0]),
		 (double)(B[1] - A[1]),
		 (double)(B[2] - A[2])};
    Vec3d AC = { (double)(C[0] - A[0]),
		 (double)(C[1] - A[1]),
		 (double)(C[2] - A[2])};

    double M_tri[21];
    mass_P2(AB, AC, M_tri);

    uint32_t pos_M = 0;
    for (uint32_t i=0; i<6; ++i)
      for (uint32_t j=0; j<=i; ++j)
	M(dof_ids[i], dof_ids[j]) += M_tri[pos_M++];
  }
}

void build_P2_stiffness_matrix(const Mesh &m,
			       const CSRPattern &P,
			       CSRMatrix &S)
{
#if USE_P2
  auto &edge2dof = m.edge2dof;
#else
  std::unordered_map<uint64_t, uint32_t> edge2dof;
#endif

  size_t vtx_count = m.vertex_count();
  size_t tri_count = m.triangle_count();
  size_t n_edges = 3 * tri_count / 2;
  assert(P.row_start.size == vtx_count + n_edges + 1);

  S.symmetric = true;
  S.rows = P.rows;
  S.cols = P.cols;
  S.nnz = P.nnz;
  S.row_start = P.row_start.data;
  S.col = P.col.data;
  S.data.resize(S.nnz);
  for (size_t i = 0; i < S.nnz; ++i) {
    S.data[i] = 0.0;
  }

  for (size_t tri=0; tri<tri_count; ++tri) {
    uint32_t v_ids[3] = {m.indices[3*tri],
			 m.indices[3*tri+1],
			 m.indices[3*tri+2]};
    std::sort(std::begin(v_ids), std::end(v_ids));
    uint32_t va = v_ids[0];
    uint32_t vb = v_ids[1];
    uint32_t vc = v_ids[2];

    uint32_t eab = edge2dof.find(pack(va, vb))->second;
    uint32_t eac = edge2dof.find(pack(va, vc))->second;
    uint32_t ebc = edge2dof.find(pack(vb, vc))->second;

    uint32_t dof_ids[6] = {va, vb, vc, eab, eac, ebc};

    Vec3f A = m.positions[va];
    Vec3f B = m.positions[vb];
    Vec3f C = m.positions[vc];

    Vec3d AB = { (double)(B[0] - A[0]),
		 (double)(B[1] - A[1]),
		 (double)(B[2] - A[2])};
    Vec3d AC = { (double)(C[0] - A[0]),
		 (double)(C[1] - A[1]),
		 (double)(C[2] - A[2])};

    double S_tri[21];
    stiffness_P2(AB, AC, S_tri);

    uint32_t pos_S = 0;
    for (uint32_t i=0; i<6; ++i)
      for (uint32_t j=0; j<=i; ++j)
	S(dof_ids[i], dof_ids[j]) += S_tri[pos_S++];
  }
}
