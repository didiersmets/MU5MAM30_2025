#include "fem_matrix.h"
#include "sparse_matrix.h"
#include "mass.h"
#include "stiffness.h"
#include "P1.h"
#include <iostream>
#include <unordered_set>
#include <vec3.h>
#include "hash_table.h"

using namespace std;

static uint64_t build_key(uint32_t i, uint32_t j)
{
    return (uint64_t(i) << 32) | uint64_t(j);
}

// -----------------------------------
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
    // dodo : clean the code, eventually write functions
    size_t nv = m.vertex_count();
    size_t ni = m.index_count();
    size_t nt = m.triangle_count();
    P.rows = nv;
    P.cols = nv;

    /* Build row_start */
    P.row_start.resize(nv + 1);
    P.row_start.set(1);
    P.row_start[0] = 0;

    for (size_t i = 0; i < ni; i++)
        P.row_start[m.indices[i] + 1] += 1;

    P.nnz = 0;
    for (size_t i = 0; i < nv; i++)
    {
        P.nnz += P.row_start[i + 1];
        P.row_start[i + 1] += P.row_start[i];
    }

    /* Build col */
    P.col.resize(P.nnz);
    P.col.set(-1);
    pair<uint32_t, uint32_t> current_pair;

    unordered_set<pair<uint32_t, uint32_t>, PairHash> last_seen;

    for (size_t i = 0; i < ni; i++)
    {
        P.col[P.row_start[m.indices[i]]] = m.indices[i];
    }

    size_t k;
    size_t i;

    for (size_t tri = 0; tri < nt; tri++)
    {
        i = 3 * tri;
        k = P.row_start[m.indices[i]] + 1;
        while (k < P.row_start[m.indices[i] + 1] && P.col[k] != -1)
            k++;

        if (k < P.row_start[m.indices[i] + 1])
        {
            /* Pair (i, i+1) */
            current_pair.first = P.row_start[m.indices[i]];
            current_pair.second = P.row_start[m.indices[i + 1]];

            if (last_seen.insert(current_pair).second)
            {
                P.col[k] = m.indices[i + 1];
                k++;
            }

            /* Pair (i, i+2) */
            current_pair.first = P.row_start[m.indices[i]];
            current_pair.second = P.row_start[m.indices[i + 2]];

            if (last_seen.insert(current_pair).second)
            {
                P.col[k] = m.indices[i + 2];
            }
        }

        k = P.row_start[m.indices[i + 1]] + 1;
        while (k < P.row_start[m.indices[i + 1] + 1] && P.col[k] != -1)
            k++;

        if (k < P.row_start[m.indices[i + 1] + 1])
        {
            /* Pair (i+1, i) */
            current_pair.first = P.row_start[m.indices[i + 1]];
            current_pair.second = P.row_start[m.indices[i]];

            if (last_seen.insert(current_pair).second)
            {
                P.col[k] = m.indices[i];
                k++;
            }

            /* Pair (i+1, i+2) */
            current_pair.first = P.row_start[m.indices[i + 1]];
            current_pair.second = P.row_start[m.indices[i + 2]];

            if (last_seen.insert(current_pair).second)
            {
                P.col[k] = m.indices[i + 2];
            }
        }

        k = P.row_start[m.indices[i + 2]] + 1;
        while (k < P.row_start[m.indices[i + 2] + 1] && P.col[k] != -1)
            k++;

        if (k < P.row_start[m.indices[i + 2] + 1])
        {
            /* Pair (i+2, i) */
            current_pair.first = P.row_start[m.indices[i + 2]];
            current_pair.second = P.row_start[m.indices[i]];

            if (last_seen.insert(current_pair).second)
            {
                P.col[k] = m.indices[i];
                k++;
            }

            /* Pair (i+2, i+2) */
            current_pair.first = P.row_start[m.indices[i + 2]];
            current_pair.second = P.row_start[m.indices[i + 1]];

            if (last_seen.insert(current_pair).second)
            {
                P.col[k] = m.indices[i + 1];
            }
        }
    }
}

// -----------------------------------
void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
    /* Data from pattern */
    M.rows = P.rows;
    M.cols = P.cols;
    M.nnz = P.nnz;

    M.row_start = static_cast<uint32_t *>(P.row_start.data);
    M.col = static_cast<uint32_t *>(P.col.data);

    M.data.resize(M.nnz);
    M.data.set(-1);

    /* Build the keys for mapping (i, j) -> k */
    HashTable<uint64_t, uint32_t> map(M.nnz);
    for (uint32_t i = 0; i < M.rows; i++)
    {
        for (size_t k = P.row_start[i]; k < P.row_start[i + 1]; k++)
        {
            uint32_t j = M.col[k];
            uint64_t key = build_key(i, j);
            map.set_at(key, k);
        }
    }

    /* Lambda map to add val in the NNZ configuration, corresponding to (i,j) */
    auto add_to_global = [&](uint32_t i, uint32_t j, double val)
    {
        uint64_t key = build_key(i, j);
        uint32_t *k = map.get(key);
        M.data[*k] = val;
    };

    double *mass_matrix = (double *)malloc(3 * 3 * sizeof(double));

    /* Build the mass matrix */
    size_t nt = m.triangle_count();
    for (size_t tri = 0; tri < nt; tri++)
    {
        uint32_t a = m.indices[3 * tri];
        uint32_t b = m.indices[3 * tri + 1];
        uint32_t c = m.indices[3 * tri + 2];

        Vec3 A = m.positions[a];
        Vec3 B = m.positions[b];
        Vec3 C = m.positions[c];

        mass(B - A, C - A, mass_matrix);

        add_to_global(a, a, mass_matrix[0]);
        add_to_global(b, b, mass_matrix[4]);
        add_to_global(c, c, mass_matrix[8]);

        add_to_global(a, b, mass_matrix[1]);
        add_to_global(b, a, mass_matrix[1]);

        add_to_global(a, c, mass_matrix[2]);
        add_to_global(c, a, mass_matrix[2]);

        add_to_global(b, c, mass_matrix[5]);
        add_to_global(c, b, mass_matrix[5]);
    }
}

// -----------------------------------
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
    /* Data from pattern */
    S.rows = P.rows;
    S.cols = P.cols;
    S.nnz = P.nnz;

    S.row_start = static_cast<uint32_t *>(P.row_start.data);
    S.col = static_cast<uint32_t *>(P.col.data);

    S.data.resize(S.nnz);
    S.data.set(-1);

    /* Build the keys for mapping (i, j) -> k */
    HashTable<uint64_t, uint32_t> map(S.nnz);
    for (uint32_t i = 0; i < S.rows; i++)
    {
        for (size_t k = P.row_start[i]; k < P.row_start[i + 1]; k++)
        {
            uint32_t j = S.col[k];
            uint64_t key = build_key(i, j);
            map.set_at(key, k);
        }
    }

    /* Lambda map to add val in the NNZ configuration, corresponding to (i,j) */
    auto add_to_global = [&](uint32_t i, uint32_t j, double val)
    {
        uint64_t key = build_key(i, j);
        uint32_t *k = map.get(key);
        S.data[*k] = val;
    };

    double *stiffness_matrix = (double *)malloc(6 * sizeof(double));

    /* Build the stiffness matrix */
    size_t nt = m.triangle_count();
    for (size_t tri = 0; tri < nt; tri++)
    {
        uint32_t a = m.indices[3 * tri];
        uint32_t b = m.indices[3 * tri + 1];
        uint32_t c = m.indices[3 * tri + 2];

        Vec3 A = m.positions[a];
        Vec3 B = m.positions[b];
        Vec3 C = m.positions[c];

        stiffness(B - A, C - A, stiffness_matrix);

        add_to_global(a, a, stiffness_matrix[0]);
        add_to_global(b, b, stiffness_matrix[1]);
        add_to_global(c, c, stiffness_matrix[2]);

        add_to_global(a, b, stiffness_matrix[3]);
        add_to_global(b, a, stiffness_matrix[3]);

        add_to_global(a, c, stiffness_matrix[5]);
        add_to_global(c, a, stiffness_matrix[5]);

        add_to_global(b, c, stiffness_matrix[4]);
        add_to_global(c, b, stiffness_matrix[4]);
    }
}