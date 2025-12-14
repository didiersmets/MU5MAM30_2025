#pragma once

#include <memory>

#include "fem_matrix.h"
#include "sparse_matrix.h"
#include "mass.h"
#include "stiffness.h"


template <typename T>
void build_P1_mass_matrix(const Mesh &m, FEMatrix &M);
template <typename T>
void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S);

template <typename T>
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P);
template <typename T>
void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M);
template <typename T>
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P,
			       CSRMatrix &S);

template <typename T>
void build_P1_SKLPattern(const Mesh &m, SKLPattern &P);
template <typename T>
void build_P1_mass_matrix(const Mesh &m, const SKLPattern &P, SKLMatrix &M);
template <typename T>
void build_P1_stiffness_matrix(const Mesh &m, const SKLPattern &P,
			       SKLMatrix &S);


template <typename T>
void build_P1_mass_matrix(const Mesh &m, FEMatrix &M) {
	double* M_loc;

	TVec3<T> AB;
	TVec3<T> AC;
	size_t tria_indices[3] = {0,0,0};  // {A_index, B_index, C_index}

	M_loc = (double *) std::malloc(9 * sizeof(double));

	for ( int tria_index = 0; tria_index < m.index_count; tria_index += 3 ) {
		tria_indices[0] = m.indices[tria_index];    // A_index
		tria_indices[1] = m.indices[tria_index+1];  // B_index
		tria_indices[2] = m.indices[tria_index+2];  // C_index

		AB = m.positions[tria_indices[1]] - m.positions[m.indices[tria_indices[0]]];
		AC = m.positions[tria_indices[2]] - m.positions[m.indices[tria_indices[0]]];

		mass(AB, AC, M_loc);

		// M[i_glob][j_glob] += M_loc[i_loc][j_loc]
		for ( int i_loc = 0; i_loc < 3; i_loc++ ) {
			for ( int j_loc = 0; j_loc < 3; j_loc++ ) {
				M[tria_indices[i_loc]][tria_indices[j_loc]] = M_loc[i_loc][j_loc];
			}
		}
	}
}

template <typename T>
void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S) {
	double* S_loc;

	TVec3<T> AB;
	TVec3<T> AC;
	size_t tria_indices[3] = {0,0,0};  // {A_index, B_index, C_index}

	M_loc = (double *) std::malloc(9 * sizeof(double));

	for ( int tria_index = 0; tria_index < m.index_count; tria_index += 3 ) {
		tria_indices[0] = m.indices[tria_index];    // A_index
		tria_indices[1] = m.indices[tria_index+1];  // B_index
		tria_indices[2] = m.indices[tria_index+2];  // C_index

		AB = m.positions[tria_indices[1]] - m.positions[m.indices[tria_indices[0]]];
		AC = m.positions[tria_indices[2]] - m.positions[m.indices[tria_indices[0]]];

		stiffness(AB, AC, S_loc);

		// M[i_glob][j_glob] += M_loc[i_loc][j_loc]
		for ( int i_loc = 0; i_loc < 3; i_loc++ ) {
			for ( int j_loc = 0; j_loc < 3; j_loc++ ) {
				S[tria_indices[i_loc]][tria_indices[j_loc]] = S_loc[i_loc][j_loc];
			}
		}
	}
}