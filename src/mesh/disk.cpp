#include <cmath>
#include <cassert>

#include "disk.h"
#include "vec3.h"

void build_disk_mesh(Mesh *m, size_t N, double R)
{
    assert(m != nullptr);

    // vertex count and triangle count as in original MyMesh implementation
    size_t vtx_count = 1 + 3 * N * (N + 1);
    size_t tri_count = 6 * N * N;

    m->positions.resize(vtx_count);
    m->indices.resize(tri_count * 3);
    m->attr.resize(vtx_count);

    // build vertices
    size_t idx = 0;
    m->positions[idx++] = Vec3(0.0f, 0.0f, 0.0f);
    for (size_t i = 1; i <= N; ++i) {
        for (size_t j = 0; j < 6 * i; ++j) {
            double angle = M_PI * j / (3.0 * i);
            double r = R * (double)i / (double)N;
            m->positions[idx++] = Vec3(static_cast<float>(r * cos(angle)),
                                       static_cast<float>(r * sin(angle)),
                                       0.0f);
        }
    }
    // ensure idx equals vtx_count
    (void)idx;

    // build triangles (ported from MyMesh implementation)
    size_t t = 0;
    // initial six triangles around center
    auto set_tri = [&](size_t a, size_t b, size_t c){
        m->indices[3 * t + 0] = static_cast<uint32_t>(a);
        m->indices[3 * t + 1] = static_cast<uint32_t>(b);
        m->indices[3 * t + 2] = static_cast<uint32_t>(c);
        ++t;
    };

    set_tri(1, 2, 0);
    set_tri(2, 3, 0);
    set_tri(3, 4, 0);
    set_tri(4, 5, 0);
    set_tri(5, 6, 0);
    set_tri(6, 1, 0);

    // further rings
    for (size_t i = 2; i <= N; ++i) {
        size_t aux = 0;
        for (size_t j = 0; j < 6 * i - 1; ++j) {
            if (j % i != 0) {
                set_tri(2 + 3 * (i - 2) * (i - 1) + aux,
                        1 + 3 * (i - 1) * (i - 2) + aux,
                        1 + 3 * (i - 1) * i + j);
                set_tri(2 + 3 * (i - 2) * (i - 1) + aux,
                        1 + 3 * (i - 1) * i + j,
                        2 + 3 * (i - 1) * i + j);
                ++aux;
            } else {
                set_tri(1 + 3 * (i - 1) * i + j,
                        2 + 3 * (i - 1) * i + j,
                        1 + 3 * (i - 2) * (i - 1) + aux);
            }
        }
        set_tri(1 + 3 * (i - 1) * i,
                3 * (i - 1) * i + 6 * i,
                2 + aux + (i - 1) * (3 * (i - 2) - 6));
        set_tri(2 + aux + (i - 1) * (3 * (i - 2) - 6),
                1 + 3 * (i - 2) * (i - 1) + aux,
                3 * (i - 1) * i + 6 * i);
    }

    // initialize attributes to zero
    for (size_t i = 0; i < vtx_count; ++i)
        m->attr[i] = 0.0f;

    // mark boundary vertices in attr (mirror MyMesh disk_boundary_edges)
    // boundary indices: 1 + 3*(N-1)*N + i, for i = 0..6*N-1
    size_t start = 1 + 3 * (N - 1) * N;
    m->boundary.resize(6 * N);
    for (size_t i = 0; i < 6 * N; ++i) {
        size_t idxb = start + i;
        if (idxb < vtx_count) {
            m->attr[idxb] = 1.0f;
            m->boundary.data[i] = static_cast<uint32_t>(idxb);
        } else {
            m->boundary.data[i] = 0;
        }
    }
}
