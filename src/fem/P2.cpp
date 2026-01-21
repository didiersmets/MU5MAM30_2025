#include <stdio.h>
#include <string.h>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass_P2.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P2.h"
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
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P)
{
    /* TOD */
}

void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
    /* TOD */
}

void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
    /* TOD */
}
