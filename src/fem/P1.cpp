#include <stdio.h>
#include <string.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "vec2.h"

using Vec2d = TVec2<double>;

static bool find_u32(uint32_t x, uint32_t *start, size_t count)
{
	for (size_t i = 0; i < count; ++i) {
		if (start[i] == x)
			return true;
	}
	return false;
}

static void stiffness2d(const Vec2d &AB, const Vec2d &AC, double *__restrict S)
{
	// compute dot products
	double ABAB = dot(AB, AB);
	double ACAC = dot(AC, AC);
	double ABAC = dot(AB, AC);

	// compute multiplier = 1 / (4 * area of the triangle)
	double mult = 0.5 / sqrt(ABAB * ACAC - ABAC * ABAC);

	// scale the dot products by the multiplier
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

static void stiffness_NS(const Vec2d &AB, const Vec2d &AC, double *__restrict S,
			 const double *u_loc, const double den, const double area)
{
	double S_start[6];

	// get the standard stiffness matrix for the triangle defined by AB and AC
	stiffness2d(AB, AC, S_start);
	
	// compute the area of the triangle
	double ABAB = dot(AB, AB);
	double ACAC = dot(AC, AC);
	double ABAC = dot(AB, AC);
	double tri_area = 0.5 * sqrt(ABAB * ACAC - ABAC * ABAC);

	// compute the dot product of grad(u) and the basis gradients: (\nabla u \cdot \nabla \phi_i)
	// compute grad(u) = sum u_loc[i] * grad(phi_i) and then compute the dot product with grad(phi_i) for i = 0,1,2
	double GraduGradPhi[3];
	GraduGradPhi[0] = u_loc[0] * S_start[0] + u_loc[1] * S_start[3] + u_loc[2] * S_start[5];
	GraduGradPhi[1] = u_loc[0] * S_start[3] + u_loc[1] * S_start[1] + u_loc[2] * S_start[4];
	GraduGradPhi[2] = u_loc[0] * S_start[5] + u_loc[1] * S_start[4] + u_loc[2] * S_start[2];

	S[0] = S_start[0] - den * den * (GraduGradPhi[0] * GraduGradPhi[0]) / tri_area;
	S[1] = S_start[1] - den * den * (GraduGradPhi[1] * GraduGradPhi[1]) / tri_area;
	S[2] = S_start[2] - den * den * (GraduGradPhi[2] * GraduGradPhi[2]) / tri_area;
	S[3] = S_start[3] - den * den * (GraduGradPhi[0] * GraduGradPhi[1]) / tri_area;
	S[4] = S_start[4] - den * den * (GraduGradPhi[1] * GraduGradPhi[2]) / tri_area;
	S[5] = S_start[5] - den * den * (GraduGradPhi[2] * GraduGradPhi[0]) / tri_area;
}

/* CSRMatrix variants */

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	size_t vtx_count = m.vertex_count(); // number of vertices
	size_t tri_count = m.triangle_count(); // number of triangles

	P.symmetric = true; 
	P.rows = P.cols = vtx_count;
	P.row_start.resize(vtx_count + 1); // row offsets --> index of the first non-zero entry in each row

	/* Upper bound on lower-triangular non-zeros (including diagonal). */
	size_t max_nnz = 3 * tri_count + vtx_count;
	P.col.resize(max_nnz); // column indices

	size_t nnz = 0;
	for (size_t a = 0; a < vtx_count; ++a) { // iterate over vertices
		P.row_start[a] = static_cast<uint32_t>(nnz);
		uint32_t *start = &P.col[nnz];
		size_t nnz_loc = 0;

		for (size_t t = 0; t < tri_count; ++t) {
			uint32_t i0 = m.indices[3 * t + 0];
			uint32_t i1 = m.indices[3 * t + 1];
			uint32_t i2 = m.indices[3 * t + 2];

			if (i0 == a) {
				if (i1 < a && !find_u32(i1, start, nnz_loc)) {
					P.col[nnz++] = i1;
					++nnz_loc;
				}
				if (i2 < a && !find_u32(i2, start, nnz_loc)) {
					P.col[nnz++] = i2;
					++nnz_loc;
				}
			} else if (i1 == a) {
				if (i2 < a && !find_u32(i2, start, nnz_loc)) {
					P.col[nnz++] = i2;
					++nnz_loc;
				}
				if (i0 < a && !find_u32(i0, start, nnz_loc)) {
					P.col[nnz++] = i0;
					++nnz_loc;
				}
			} else if (i2 == a) {
				if (i0 < a && !find_u32(i0, start, nnz_loc)) {
					P.col[nnz++] = i0;
					++nnz_loc;
				}
				if (i1 < a && !find_u32(i1, start, nnz_loc)) {
					P.col[nnz++] = i1;
					++nnz_loc;
				}
			}
		}

		P.col[nnz++] = static_cast<uint32_t>(a);
	}

	P.row_start[vtx_count] = static_cast<uint32_t>(nnz);
	P.col.resize(nnz);
	P.col.shrink_to_fit();
	P.nnz = nnz;

	for (size_t a = 0; a < vtx_count; ++a) {
		uint32_t *to_sort = &P.col[P.row_start[a]];
		size_t count = P.row_start[a + 1] - P.row_start[a];
		for (size_t k = 1; k < count; ++k) {
			size_t j = k - 1;
			while (j && (to_sort[j] > to_sort[j + 1])) {
				uint32_t tmp = to_sort[j];
				to_sort[j] = to_sort[j + 1];
				to_sort[j + 1] = tmp;
				--j;
			}
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
}

static void init_CSRMatrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	S.symmetric = true;
	S.rows = S.cols = m.vertex_count();
	S.nnz = P.col.size;
	S.row_start = P.row_start.data;
	S.col = P.col.data;
	S.data.resize(S.nnz);
	for (size_t i = 0; i < S.nnz; ++i)
		S.data[i] = 0.0;
}


void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	assert(P.row_start.size == m.vertex_count() + 1);
	init_CSRMatrix(m, P, S);

	for (size_t t = 0; t < m.triangle_count(); ++t) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		Vec2d A = { static_cast<double>(m.positions[a].x), static_cast<double>(m.positions[a].y) };
		Vec2d B = { static_cast<double>(m.positions[b].x), static_cast<double>(m.positions[b].y) };
		Vec2d C = { static_cast<double>(m.positions[c].x), static_cast<double>(m.positions[c].y) };
		Vec2d AB = { B[0] - A[0], B[1] - A[1] };
		Vec2d AC = { C[0] - A[0], C[1] - A[1] };
		double Sloc[6];
		stiffness2d(AB, AC, Sloc);
		S(a, a) += Sloc[0];
		S(b, b) += Sloc[1];
		S(c, c) += Sloc[2];
		S(a > b ? a : b, a > b ? b : a) += Sloc[3];
		S(b > c ? b : c, b > c ? c : b) += Sloc[4];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5];
	}
}


void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S,
				const double *den)
{
	assert(P.row_start.size == m.vertex_count() + 1);
	init_CSRMatrix(m, P, S);

	for (size_t t = 0; t < m.triangle_count(); ++t) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		Vec2d A = { static_cast<double>(m.positions[a].x), static_cast<double>(m.positions[a].y) };
		Vec2d B = { static_cast<double>(m.positions[b].x), static_cast<double>(m.positions[b].y) };
		Vec2d C = { static_cast<double>(m.positions[c].x), static_cast<double>(m.positions[c].y) };
		Vec2d AB = { B[0] - A[0], B[1] - A[1] };
		Vec2d AC = { C[0] - A[0], C[1] - A[1] };
		double Sloc[6];
		stiffness2d(AB, AC, Sloc);
		double w = den[t];
		S(a, a) += w * Sloc[0];
		S(b, b) += w * Sloc[1];
		S(c, c) += w * Sloc[2];
		S(a > b ? a : b, a > b ? b : a) += w * Sloc[3];
		S(b > c ? b : c, b > c ? c : b) += w * Sloc[4];
		S(c > a ? c : a, c > a ? a : c) += w * Sloc[5];
	}
}

void build_P1_stiffness_matrix_NS(const Mesh &m, const CSRPattern &P,
				  CSRMatrix &S, const double *den,
				  const double *u,
				  double area
				)
{
	/*
	implements the computation of global stiffness matrix J(u)
	S is global stiffness matrix only in CSR format
	P is the pattern of the stiffness matrix S, built from the mesh m
	m is the mesh of the problem, containing the vertices and triangles
	u is the current solution, used to compute the local stiffness matrix for each triangle
	den is the denominator of the energy functional, used to compute the local stiffness matrix for each triangle
	*/
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

	// for each triangle we compute local stiffness matrix and add it to the global stiffness matrix
	for (size_t t = 0; t < tri_count; ++t) {

		// get vertices of the triangle
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];

		// get coordinates of the vertices and compute AB and AC
		Vec2d A = { static_cast<double>(m.positions[a].x),
			    static_cast<double>(m.positions[a].y) };
		Vec2d B = { static_cast<double>(m.positions[b].x),
			    static_cast<double>(m.positions[b].y) };
		Vec2d C = { static_cast<double>(m.positions[c].x),
			    static_cast<double>(m.positions[c].y) };
		Vec2d AB = { B[0] - A[0], B[1] - A[1] };
		Vec2d AC = { C[0] - A[0], C[1] - A[1] };

		// get local solution at the vertices of the triangle
		double u_loc[3] = { u[a], u[b], u[c] };
		double Sloc[6];

		// compute local stiffnes matrix
		stiffness_NS(AB, AC, Sloc, u_loc, den[t], area);

		// add local stiffness matrix to global stiffness matrix based on the pattern P
		// diagonal entries
		S(a, a) += Sloc[0] * den[t];
		S(b, b) += Sloc[1] * den[t];
		S(c, c) += Sloc[2] * den[t];

		// off-diagonal entries
		// since CSR is upper triangular we only add the entries for (a, b), (b, c) and (c, a) if they are in the pattern P
		S(a > b ? a : b, a > b ? b : a) += Sloc[3] * den[t]; 
		S(b > c ? b : c, b > c ? c : b) += Sloc[4] * den[t];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5] * den[t];
	}
}

void build_P1_rhs_NS(const Mesh &m, const double *den, const double *u,
		     TArray<double> &rhs)
{
	memset(rhs.data, 0, m.vertex_count() * sizeof(double));

	// loop over all triangles
	for (size_t t = 0; t < m.triangle_count(); ++t) {

		// get vertices of the triangle
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];

		// get coordinates of the vertices and compute AB and AC
		Vec2d A = { static_cast<double>(m.positions[a].x),
			    static_cast<double>(m.positions[a].y) };
		Vec2d B = { static_cast<double>(m.positions[b].x),
			    static_cast<double>(m.positions[b].y) };
		Vec2d C = { static_cast<double>(m.positions[c].x),
			    static_cast<double>(m.positions[c].y) };
		Vec2d AB = { B[0] - A[0], B[1] - A[1] };
		Vec2d AC = { C[0] - A[0], C[1] - A[1] };

		double S_loc[6];

		stiffness2d(AB, AC, S_loc);

		// scale by the denominator of the nonlinearity
		rhs[a] -= den[t] * (u[a] * S_loc[0] + u[b] * S_loc[3] + u[c] * S_loc[5]);
		rhs[b] -= den[t] * (u[a] * S_loc[3] + u[b] * S_loc[1] + u[c] * S_loc[4]);
		rhs[c] -= den[t] * (u[a] * S_loc[5] + u[b] * S_loc[4] + u[c] * S_loc[2]);
	}
}

/* FEMatrix variants */
/*
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
	*/

/*
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
*/

