#pragma once

#include "fem_matrix.h"
#include "sparse_matrix.h"

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
	
}