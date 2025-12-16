#pragma once

#include "fem_matrix.h"
#include "sparse_matrix.h"

/* CSR variants */
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P);
static uint64_t build_key(uint32_t i, uint32_t j);

void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S);

void build_P1_SKLPattern(const Mesh &m, SKLPattern &P);
void build_P1_mass_matrix(const Mesh &m, const SKLPattern &P, SKLMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, const SKLPattern &P,
							   SKLMatrix &S);

struct PairHash
{
	std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const noexcept
	{
		return std::hash<uint32_t>{}(p.first) ^ (std::hash<uint32_t>{}(p.second) << 1);
	}
};
