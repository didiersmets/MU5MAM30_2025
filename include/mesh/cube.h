#pragma once
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"

static void load_cube_vertices(Vec3 *pos, size_t subdiv)

{
    size_t n= subdiv +1;

    size_t voff[6]
    for (int f; f<6; f++) {
        voff[f]=f*POW2(n);
    }

    float *dir = (float *)safe_malloc(n * sizeof(float));
    float *rev = (float *)safe_malloc(n * sizeof(float));
    for (size_t i = 0; i < n; ++i) {
        dir[i] = (2 * (float)i - subdiv) / subdiv;
    }

    for (size_t i = 0; i < n; ++i){
        rev[i] = dir[subdiv-i];
    }
    for (size_t y = 0; y < n; ++y) {
        for (size_t x = 0; x < n; ++x) {
            pos[voff[0]++] = { dir[x], -1, dir[y] }; /* Front  */
            pos[voff[1]++] = { rev[x], +1, dir[y] }; /* Back   */
            pos[voff[2]++] = { -1, rev[x], dir[y] }; /* Left   */
            pos[voff[3]++] = { +1, dir[x], dir[y] }; /* Right  */
            pos[voff[4]++] = { dir[x], rev[y], -1 }; /* Bottom */
            pos[voff[5]++] = { dir[x], dir[y], +1 }; /* Top    */
        }
    }
    free(dir);
    free(rev);

}

static void load_cube_indices(uint32_t *idx, size_t subdiv) {
    size_t n= subdiv +1;
    for (size_t f = 0; f<6; f++) {
        size_t offset = f*POW2(n);
        for (size_t i=0; i<subdiv;i++) {
            for (size_t j=0; j<subdiv;j++) {
                int base = i*n+j+offset;
                *idx++ = base;
                *idx++ = base+1;
                *idx++ = base +1+n;

                *idx++ = base;
                *idx++ = base +1+n;
                *idx++ = base +n;

            }
        }
    }
}

int load_cube(Mesh &m, size_t subdiv) {

    size_t n = subdiv + 1;

    /* Reserve vertices and indices */
    m.positions.resize(6 * POW2(n));
    m.indices.resize(36 * POW2(subdiv));

    /* First build vertices as six unattached faces of n^2 vertices each */
    load_cube_vertices(m.positions.data, subdiv);

    /* Build corresponding triangulation indices */
    load_cube_indices(m.indices.data, subdiv);

    /* Finally attach faces between themselves */
    remove_duplicate_vertices(m);

    return (0);
}
