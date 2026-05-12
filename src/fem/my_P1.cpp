#include <stdio.h>
#include <string.h>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

#include "my_P1.h"
#include "my_adjacency.h"
#include "my_mesh.h"
#include "my_sparse_matrix.h"
#include "vec2.h"

using Vec2d = TVec2<double>;

static bool find(uint32_t x, uint32_t *start, size_t count)
{
	for (size_t i = 0; i < count; ++i) {
		if (start[i] == x)
			return true;
	}
	return false;
}

void build_P1_CSRPattern(const MyMesh &m, MyCSRPattern &P)
{
	size_t vtx_count = m.vtx_count;
	size_t tri_count = m.tri_count;

	P.row_start.resize(vtx_count + 1);

	MyVTAdjacency adj(m);

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
		uint32_t *to_sort = &P.col[P.row_start[a]];
		size_t count = P.row_start[a + 1] - P.row_start[a];
		/* Insertion sort */
		for (size_t k = 1; k < count; ++k) {
			size_t j = k - 1;
			while (j && (to_sort[j] > to_sort[j + 1])) {
				uint32_t tmp = to_sort[j];
				to_sort[j] = to_sort[j + 1];
				to_sort[j + 1] = tmp;
				j--;
			}
		}
	}
}

static void stiffness(const Vec2d &AB, const Vec2d &AC, double *__restrict S);
static void stiffness_NS(const Vec2d &AB, const Vec2d &AC, double *__restrict S, const double *u, const double den);

void inline stiffness(const Vec2d &AB, const Vec2d &AC, double *__restrict S)
{
	/*
	compute the local stiffness matrix for the triangle defined by the vectors AB and AC
	only compute triangular matrix of the 3x3 and s
	*/

	// compute the dot products needed for the stiffness matrix
	double ABAB = dot(AB,AB);
	double ACAC = dot(AC,AC);
	double ABAC = dot(AB, AC);

	// multiplier is 1/4 of the area of the triangle
	double mult = 0.5 / sqrt(ABAB * ACAC - ABAC * ABAC); 

	// scale the dot products by the multiplier to get the entries of the stiffness matrix
	ABAB *= mult;
	ACAC *= mult;
	ABAC *= mult;

	// the entries of the local stiffness matrix S are computed from the scaled dot products
	// S_ij = T_area * (grad phi_i . grad phi_j) = T_area * (1/(2*T_area) * (edge opposite to i) . (edge opposite to j)) = 1/4 * (edge opposite to i) . (edge opposite to j)
	
	S[0] = ACAC - 2 * ABAC + ABAB;  //S(A,A)= dot(grad phi_0, grad phi_0) = 1/(4*T_area) * (edge BC . edge BC) = 1/(4*T_area) * (AC - AB) . (AC - AB) = 1/(4*T_area) * (AB.AB + AC.AC - 2*AB.AC)
	S[1] = ACAC; 					// S(B,B)= dot(grad phi_1, grad phi_1) = 1/(4*T_area) * (edge CA . edge CA) = 1/(4*T_area) * -AC . -AC
	S[2] = ABAB; 					// S(C,C)= dot(grad phi_2, grad phi_2) = 1/(4*T_area) * (edge AB . edge AB) = 1/(4*T_area) * AB . AB
	S[3] = ABAC - ACAC; 			// S(A,B)= dot(grad phi_0, grad phi_1) = 1/(4*T_area) * (edge BC . edge CA) = 1/(4*T_area) * (AC - AB) . AC = 1/(4*T_area) * (AC.AC - AB.AC)
	/* Note the chosen order : (B,C)-> 4 and (C,A) -> 5 */
	S[4] = -ABAC; 					// S(B,C)= dot(grad phi_1, grad phi_2) = 1/(4*T_area) * (edge CA . edge AB) = 1/(4*T_area) * -AC . AB = -1/(4*T_area) * AB.AC
	S[5] = ABAC - ABAB; 			// S(C,A)= dot(grad phi_2, grad phi_0) = 1/(4*T_area) * (edge AB . edge BC) = 1/(4*T_area) * AB . (AC - AB) = 1/(4*T_area) * (AB.AB - AB.AC)

}

static void stiffness_NS(const Vec2d &AB, const Vec2d &AC, double *__restrict S, const double *u_loc, const double den){
	// compute the standard linear stiffness matrix
	// S_start = \int_T grad(phi_i) dot grad(phi_j)
	double S_start[6];
	stiffness(AB, AC, S_start); 
	
	// compute the dot product of grad(u) and the basis gradients: (\nabla u \cdot \nabla \phi_i)
	// compute grad(u) = sum u_loc[i] * grad(phi_i) and then compute the dot product with grad(phi_i) for i = 0,1,2
	double GraduGradPhi[3];
	GraduGradPhi[0] = u_loc[0]*S_start[0] + u_loc[1]*S_start[3] + u_loc[2]*S_start[5];
	GraduGradPhi[1] = u_loc[0]*S_start[3] + u_loc[1]*S_start[1] + u_loc[2]*S_start[4];
	GraduGradPhi[2] = u_loc[0]*S_start[5] + u_loc[1]*S_start[4] + u_loc[2]*S_start[2];

	// the jacobian stiffness matrix becomes S_ij = grad(phi_i) dot grad(phi_j) - den**2 * (grad(u) dot grad(phi_i)) * (grad(u) dot grad(phi_j))
	// den is the denominator of the energy functional, used to scale the non-linear part of the stiffness matrix
	S[0] = S_start[0] - den*den*(GraduGradPhi[0]*GraduGradPhi[0]);
	S[1] = S_start[1] - den*den*(GraduGradPhi[1]*GraduGradPhi[1]);
	S[2] = S_start[2] - den*den*(GraduGradPhi[2]*GraduGradPhi[2]);
	S[3] = S_start[3] - den*den*(GraduGradPhi[0]*GraduGradPhi[1]);
	S[4] = S_start[4] - den*den*(GraduGradPhi[1]*GraduGradPhi[2]);
	S[5] = S_start[5] - den*den*(GraduGradPhi[2]*GraduGradPhi[0]);
}

void build_P1_stiffness_matrix(const MyMesh &m, const MyCSRPattern &P, MyCSRMatrix &S, bool modified, const double *den) 
{

	/*
	implementation of the computation of global stiffness matrix S for the functional
	*/
	size_t vtx_count = m.vtx_count;
	size_t tri_count = m.tri_count;
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
	size_t a,b,c;
	Vec2d A, B, C, AB, AC;
	double Sloc[6];
	for (size_t t = 0; t < tri_count; ++t) {
		a = m.triangles[t].x;
		b = m.triangles[t].y;
		c = m.triangles[t].z;

		A = m.vertices[a];
		B = m.vertices[b];
		C = m.vertices[c];
		AB = {(double)B[0] - (double)A[0],
			    (double)B[1] - (double)A[1]};
		AC = {(double)C[0] - (double)A[0],
			    (double)C[1] - (double)A[1]};
		
		stiffness(AB, AC, Sloc);
		if (modified) {
			Sloc[0] *= den[t];
			Sloc[1] *= den[t];
			Sloc[2] *= den[t];
			Sloc[3] *= den[t];
			Sloc[4] *= den[t];
			Sloc[5] *= den[t];
		}

		S(a, a) += Sloc[0];
		S(b, b) += Sloc[1];
		S(c, c) += Sloc[2];
		S(a > b ? a : b, a > b ? b : a) += Sloc[3];
		S(b > c ? b : c, b > c ? c : b) += Sloc[4];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5];
	}
}



void build_P1_stiffness_matrix_NS(const MyMesh &m, const MyCSRPattern &P, MyCSRMatrix &S, const double *den, const double *u) 
{
	/*
	implements the computation of global stiffness matrix J(u)
	S is global stiffness matrix only in CSR format
	P is the pattern of the stiffness matrix S, built from the mesh m
	m is the mesh of the problem, containing the vertices and triangles
	u is the current solution, used to compute the local stiffness matrix for each triangle
	den is the denominator of the energy functional, used to compute the local stiffness matrix for each triangle
	*/
	size_t vtx_count = m.vtx_count;
	size_t tri_count = m.tri_count;
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

	size_t a,b,c;
	Vec2d A, B, C, AB, AC;
	double Sloc[6], u_loc[3];

	// for each triangle compute the local stiffness matrix and add it to the global one
	for (size_t t = 0; t < tri_count; ++t) {

		// get the indices of the vertices of the triangle
		a = m.triangles[t].x;
		b = m.triangles[t].y;
		c = m.triangles[t].z;

		// get the coordinates of the vertices of the triangle
		A = m.vertices[a];
		B = m.vertices[b];
		C = m.vertices[c];

		// compute the vectors AB and AC
		AB = {(double)B[0] - (double)A[0],
			    (double)B[1] - (double)A[1]};
		AC = {(double)C[0] - (double)A[0],
			    (double)C[1] - (double)A[1]};
	
		// get the values of u at the vertices of the triangle
		u_loc[0] = u[a];
		u_loc[1] = u[b];
		u_loc[2] = u[c];

		// compute local stiffness matrix Sloc
		stiffness_NS(AB, AC, Sloc, u_loc, den[t]);
		
		// add Sloc to the global stiffness matrix S scaled by den[t]
		S(a, a) += Sloc[0]*den[t];
		S(b, b) += Sloc[1]*den[t];
		S(c, c) += Sloc[2]*den[t];

		// Add off-diagonal elements to the global symmetric matrix S
		// We use max/min to ensure we only write to the lower triangular part.
		// we only store S in CSR format with the convention that S(i,j) with i >= j is stored, so we need to check the order of the indices
		S(a > b ? a : b, a > b ? b : a) += Sloc[3]*den[t];
		S(b > c ? b : c, b > c ? c : b) += Sloc[4]*den[t];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5]*den[t];
	}
}

void build_P1_rhs_NS(const MyMesh &m, const double *den, const double *u, TArray<double> &rhs){

	/*
	build right hand side for the non linear problem, which is the gradient of the energy functional
	we compute the local contribution of each triangle and add it to the global right hand side
	*/
	
	memset(rhs.data, 0, m.vtx_count * sizeof(double));
	double S_loc[6];
	size_t a,b,c;
	Vec2d A, B, C, AB, AC;
	
	// loop over all triangles
	for(size_t t = 0; t < m.tri_count; t++){

		// get global indices of the vertices of the triangle
		a = m.triangles[t].x;
		b = m.triangles[t].y;
		c = m.triangles[t].z;
		
		// get the coordinates of the vertices of the triangle
		A = m.vertices[a];
		B = m.vertices[b];
		C = m.vertices[c];

		// compute the vectors AB and AC in x and y
		AB = {(double)B[0] - (double)A[0],
			    (double)B[1] - (double)A[1]};
		AC = {(double)C[0] - (double)A[0],
			    (double)C[1] - (double)A[1]};

		// compute local stiffness matrix for the triangle, used to compute the contribution of the triangle to the right hand side
		stiffness(AB, AC, S_loc);

		// ad contribution of trinagle to right hand side
		rhs[a] -= den[t]*(u[a]*S_loc[0] +u[b]*S_loc[3] + u[c]*S_loc[5]); // row a = den * (S(a,a)*u[a] + S(a,b)*u[b] + S(a,c)*u[c])
		rhs[b] -= den[t]*(u[a]*S_loc[3] +u[b]*S_loc[1] + u[c]*S_loc[4]); // row b = den * (S(b,a)*u[a] + S(b,b)*u[b] + S(b,c)*u[c]) but S(b,a) = S(a,b) = S_loc[3]
		rhs[c] -= den[t]*(u[a]*S_loc[5] +u[b]*S_loc[4] + u[c]*S_loc[2]); // row c = den * (S(c,a)*u[a] + S(c,b)*u[b] + S(c,c)*u[c]) but S(c,a) = S(a,c) = S_loc[5]
	}

}