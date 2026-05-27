#include <cmath>
#include <cassert>

#include "square.h"
#include "vec3.h"

void build_square_mesh(Mesh *m, size_t N, double a, double b)
{
    assert(m != nullptr);
    size_t V = N + 1;
    size_t vtx_count = V * V;
    size_t tri_count = 2 * N * N;

    m->positions.resize(vtx_count);
    m->indices.resize(tri_count * 3);
    m->attr.resize(vtx_count);

    // vertices
    for (size_t i = 0; i < V; ++i) {
        for (size_t j = 0; j < V; ++j) {
            size_t idx = V * i + j;
            m->positions[idx] = Vec3(static_cast<float>(a * i / (double)N),
                                     static_cast<float>(b * j / (double)N),
                                     0.0f);
            m->attr[idx] = 0.0f;
        }
    }

    // triangles
    size_t t = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            size_t v = V * i + j;
            // triangle 1: v, v+1, v+1+V
            m->indices[3 * (t) + 0] = static_cast<uint32_t>(v);
            m->indices[3 * (t) + 1] = static_cast<uint32_t>(v + 1);
            m->indices[3 * (t) + 2] = static_cast<uint32_t>(v + 1 + V);
            ++t;
            // triangle 2: v, v+1+V, v+V
            m->indices[3 * (t) + 0] = static_cast<uint32_t>(v);
            m->indices[3 * (t) + 1] = static_cast<uint32_t>(v + 1 + V);
            m->indices[3 * (t) + 2] = static_cast<uint32_t>(v + V);
            ++t;
        }
    }

    // mark boundary vertices in attr (mirror MyMesh boundary array)
    // boundary order in MyMesh: 0, N, V*N, V*N+N, then for i=1..N-1: i, V*i, V*i+N, V*N+i
    m->attr[0] = 1.0f;
    m->attr[N] = 1.0f;
    m->attr[V * N] = 1.0f;
    m->attr[V * N + N] = 1.0f;
    m->boundary.resize(4 * N);
    size_t k = 0;
    m->boundary.data[k++] = 0;
    m->boundary.data[k++] = static_cast<uint32_t>(N);
    m->boundary.data[k++] = static_cast<uint32_t>(V * N);
    m->boundary.data[k++] = static_cast<uint32_t>(V * N + N);
    for (size_t i = 1; i < N; ++i) {
        m->attr[i] = 1.0f;               // bottom side
        m->attr[V * i] = 1.0f;           // right side
        m->attr[V * i + N] = 1.0f;       // left side
        m->attr[V * N + i] = 1.0f;       // top side
        m->boundary.data[k++] = static_cast<uint32_t>(i);
        m->boundary.data[k++] = static_cast<uint32_t>(V * i);
        m->boundary.data[k++] = static_cast<uint32_t>(V * i + N);
        m->boundary.data[k++] = static_cast<uint32_t>(V * N + i);
    }
}
