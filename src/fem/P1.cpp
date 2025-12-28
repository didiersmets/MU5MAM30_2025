#include "P1.h"
#include "mass.h"
#include "stiffness.h"

void qsort(TArray<uint32_t> &a, size_t start, size_t end){
	if (start >= end-1) return;
	// todo implement quicksort
	uint32_t pivot = a[start];
	size_t mid = start+1;
	for (size_t i=start+1;i<end;i++){
		if (a[i]<=pivot){
			uint32_t temp = a[i];
			a[i] = a[mid];
			a[mid] = temp;
			temp++;
		}
	}
	qsort(a,start,mid);
	qsort(a,mid,end);
}

void build_P1_CSRPattern(const Mesh &m, CSRPattern &pattern){
	// TODO : rework : ugly
	pattern.rows = m.vertex_count();
	pattern.cols = pattern.rows;
	TArray<uint32_t> row_start (pattern.rows+1,1);
	for(size_t t=0;t<m.index_count();t+=3){
		for (int e=0;e<3;e++){
			uint32_t v = m.indices[t+e];
			row_start[v+1]++;
		}
	}
	// cumsum to get row_start
	for (size_t i=0;i<pattern.rows;i++){
		row_start[i+1] += row_start[i];
	}
	pattern.nnz = row_start[pattern.rows];

	TArray<uint32_t> col(pattern.nnz);
	TArray<uint32_t> counter(pattern.rows, 1);

	//diagonal
	for (size_t i = 0;i<pattern.rows;i++) col[row_start[i]]=i;

	for (size_t t=0;t<m.index_count();t+=3){
		for (int e=0;e<3;e++){
			uint32_t v1 = m.indices[t+e];
			uint32_t v2 = m.indices[t+(e+1)%3];
			uint32_t start1 = row_start[v1];
			bool found = false;
			for(size_t i=0;i<counter[v1];i++){
				if (col[start1 + i]==v2){
					found = true;
				}
			}
			if (!found){
				col[start1 + counter[v1]] = v2;
				counter[v1]++;
			}
			uint32_t start2 = row_start[v2];
			found = false;
			for(size_t i=0;i<counter[v2];i++){
				if (col[start2 + i]==v1){
					found = true;
				}
			}
			if (!found){
				col[start2 + counter[v2]] = v1;
				counter[v2]++;
			}
		}
	}
	for(size_t i = 0;i<pattern.rows;i++) qsort(col,row_start[i],row_start[i+1]);
	pattern.row_start = std::move(row_start);
	pattern.col = std::move(col);
}

void build_P1_mass_matrix(const Mesh &m, CSRMatrix &M){
	for (size_t t=0;t<m.index_count();t+=3){
		uint32_t Ai = m.indices[t];
		uint32_t Bi = m.indices[t+1];
		uint32_t Ci = m.indices[t+2];
		Vec3 A = m.positions[Ai];
		Vec3 B = m.positions[Bi];
		Vec3 C = m.positions[Ci];
		Vec3 AB = B-A;
		Vec3 AC = C-A;
		MassCoef coef = mass(AB,AC);
		M(Ai,Ai)+=coef.diag;
		M(Bi,Bi)+=coef.diag;
		M(Ci,Ci)+=coef.diag;
		M(Ai,Bi)+=coef.offdiag;
		M(Bi,Ci)+=coef.offdiag;
		M(Ci,Ai)+=coef.offdiag;
		M(Bi,Ai)+=coef.offdiag;
		M(Ci,Bi)+=coef.offdiag;
		M(Ai,Ci)+=coef.offdiag;
	}
}
void build_P1_stiffness_matrix(const Mesh &m, CSRMatrix &S){
	for (size_t t=0;t<m.index_count();t+=3){
		uint32_t Ai = m.indices[t];
		uint32_t Bi = m.indices[t+1];
		uint32_t Ci = m.indices[t+2];
		Vec3 A = m.positions[Ai];
		Vec3 B = m.positions[Bi];
		Vec3 C = m.positions[Ci];
		Vec3 AB = B-A;
		Vec3 AC = C-A;
		StiffnessCoef coef = stiffness(AB,AC);
		S(Ai,Ai)+=coef.S00;
		S(Bi,Bi)+=coef.S11;
		S(Ci,Ci)+=coef.S22;
		S(Ai,Bi)+=coef.S01;
		S(Bi,Ci)+=coef.S12;
		S(Ci,Ai)+=coef.S20;
		S(Bi,Ai)+=coef.S01;
		S(Ci,Bi)+=coef.S12;
		S(Ai,Ci)+=coef.S20;
	}
}
