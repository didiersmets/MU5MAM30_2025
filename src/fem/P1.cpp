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

bool find(const TArray<uint32_t> & arr, uint32_t target,  size_t start, size_t stop){
	assert(start <= stop && stop <= sizeof(arr));

	for(size_t i = start; i < stop; i ++){
		if(arr[i] == target){return true;}
	}
	return false;
}

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	P.symmetric = true;

	size_t num_vtx = m.vertex_count();
	P.rows = num_vtx;
	P.cols = num_vtx;
	
	VTAdjacency adj(m);

	size_t nnz_max  = 3 * m.triangle_count() + num_vtx;

	//Initilalize 
	size_t nnz = 0;
	P.row_start[0] = 0;
	P.col.resize(nnz_max);

	for(size_t i = 0; i < num_vtx; i ++ ){
		uint32_t start = adj.offset[i];
		uint32_t stop = start + adj.degree[i];
		
		for(uint32_t j = start; j < stop; j++){
			//look who is adjacent to current index;
			
			uint32_t a = adj.vtri[j].prev;
			uint32_t b = adj.vtri[j].next;

			size_t row_start = P.row_start[i];

			//see if we already added them to the matrix
			if(a < i && !find(P.col, a, row_start, nnz)){
				P.col[nnz++] = a;
			}
			if(b < i && !find(P.col, b, row_start, nnz)){
				P.col[nnz++] = b;
			}	
		}
		P.col[nnz++] = i;  // set the diagonal, this will always be the last one anyway.
		P.row_start[i+1] = nnz;
	}
	P.nnz = nnz;
	P.col.resize(P.nnz);
	P.col.shrink_to_fit();

	//reorder the indicies
	uint32_t stop = 0;
	for(size_t i = 0; i < num_vtx; i ++){
		uint32_t start = stop;
		stop = P.row_start[i+1];

		for(size_t j = start+1; j < stop; j++){
			for(size_t k = start; k <j; k++){
				if(P.col[j] < P.col[k]){
					uint32_t tmp = P.col[j];
					P.col[j] = P.col[k];
					P.col[k] = tmp;
				}
			}
		}
	}
}

void add_to_mass_matrix(uint32_t start, uint32_t stop, uint32_t ind1, uint32_t ind2, uint32_t ind3, 
					uint32_t *col, TArray<double>& data, double* M ){

// funcion to add att the right place in the mass matrix

	if(ind2 < ind1){
		for(uint32_t j = start; j < stop; j ++){
			if(col[j] == ind2){
				data[j] += M[1];
			}
		}
		}
	if(ind3 < ind1){
		for(uint32_t j = start; j < stop; j ++){
			if(col[j] == ind3){
				data[j] += M[1];
			}
		}
		}
		data[stop-1] += M[0];	
}



void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
	size_t vtx_count = m.vertex_count();
	//size_t tri_count = m.triangle_count();
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
	


	for(size_t i = 0; i < m.index_count(); i+=3){
		//indicies of current triangle
		uint32_t ind1 = m.indices[i];
		uint32_t ind2 = m.indices[i+1];
		uint32_t ind3 = m.indices[i+2];
		//position of current vectors
		Vec3 v1 = m.positions[ind1];
		Vec3 v2 = m.positions[ind2];
		Vec3 v3 = m.positions[ind3];
		
		//Compute local mass matrix
		double mass_mat[2];
		Vec3d AB = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
		Vec3d AC = {v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]};
		mass(AB, AC, mass_mat);
		uint32_t start, stop;

		//start of the row v1
		 start = P.row_start[ind1];
		 stop = P.row_start[ind1 +1];

		add_to_mass_matrix( start,  stop,  ind1,  ind2,  ind3, M.col, M.data, mass_mat );

		//start of the row v1
		 start = P.row_start[ind2];
		 stop = P.row_start[ind2 +1];

		add_to_mass_matrix( start,  stop,  ind2,  ind1,  ind3, M.col, M.data, mass_mat );
		
		//start of the row v1
		start = P.row_start[ind3];
		stop = P.row_start[ind3 +1];

		add_to_mass_matrix( start,  stop,  ind3,  ind2,  ind1, M.col, M.data, mass_mat );
	
	}
}

void add_to_stiffness_matrix(uint32_t start, uint32_t stop, uint32_t ind1, uint32_t ind2, uint32_t ind3, 
					uint32_t *col, TArray<double>& data, double first_off_diag, double second_off_diag, double diag ){

// funcion to add att the right place in the stifness matrix

	if(ind2 < ind1){
		for(uint32_t j = start; j < stop; j ++){
			if(col[j] == ind2){
				data[j] += first_off_diag;
			}
		}
		}
	if(ind3 < ind1){
		for(uint32_t j = start; j < stop; j ++){
			if(col[j] == ind3){
				data[j] += second_off_diag;
			}
		}
		}
		data[stop-1] += diag;	
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	size_t vtx_count = m.vertex_count();
	//size_t tri_count = m.triangle_count();
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

		for(size_t i = 0; i < m.index_count(); i+=3){
		//indicies of current triangle
		uint32_t ind1 = m.indices[i];
		uint32_t ind2 = m.indices[i+1];
		uint32_t ind3 = m.indices[i+2];
		//position of current vectors
		Vec3 v1 = m.positions[ind1];
		Vec3 v2 = m.positions[ind2];
		Vec3 v3 = m.positions[ind3];
		
		//Compute local stiffness matrix
		double stiff_mat[6];
		Vec3d AB = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
		Vec3d AC = {v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]};
		stiffness(AB, AC, stiff_mat);
		uint32_t start, stop;

		//start of the row v1
		 start = P.row_start[ind1];
		 stop = P.row_start[ind1 +1];

		add_to_stiffness_matrix( start,  stop,  ind1,  ind2,  ind3, S.col, S.data,
			 stiff_mat[3], stiff_mat[4], stiff_mat[0] );

		//start of the row v2
		 start = P.row_start[ind2];
		 stop = P.row_start[ind2 +1];

		add_to_stiffness_matrix( start,  stop,  ind2,  ind1,  ind3, S.col, S.data,
			 stiff_mat[3], stiff_mat[5], stiff_mat[1] );
		
			 //start of the row v3
		start = P.row_start[ind3];
		stop = P.row_start[ind3 +1];

		add_to_stiffness_matrix( start,  stop,  ind3,  ind1,  ind2, S.col, S.data,
			 stiff_mat[4], stiff_mat[5], stiff_mat[2] );	
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
