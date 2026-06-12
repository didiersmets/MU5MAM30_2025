#include <cmath>
#include <cassert>

#include "square.h"
#include "vec3.h"

void build_square_mesh(Mesh *m, size_t N, double a, double b)
{
    assert(m != nullptr);
    size_t V = N + 1; // number of vertices per side
    size_t vtx_count = V * V; // total number of vertices
    size_t tri_count = 2 * N * N; // total number of triangles (2 per square cell)

    m->positions.resize(vtx_count);
    m->indices.resize(tri_count * 3); // 3 indices per triangle
    m->attr.resize(vtx_count);

    // vertices
    for (size_t i = 0; i < V; ++i) { // iterate over rows
        for (size_t j = 0; j < V; ++j) { // iterate over columns
            size_t idx = V * i + j; // vertex index in row-major order
            m->positions[idx] = Vec3(static_cast<float>(a * i / (double)N), // x coordinate, a is the width of the square
                                     static_cast<float>(b * j / (double)N), // y coordinate, b is the height of the square
                                     0.0f);
            m->attr[idx] = 0.0f;
        }
    }

    // triangles
    size_t t = 0;
    for (size_t i = 0; i < N; ++i) { // iterate over rows of cells
        for (size_t j = 0; j < N; ++j) { // iterate over columns of cells
            size_t v = V * i + j; // row major index of the bottom-left vertex of the current cell
            // triangle 1: v, v+1, v+1+V, 
            m->indices[3 * (t) + 0] = static_cast<uint32_t>(v); // bottom-left vertex
            m->indices[3 * (t) + 1] = static_cast<uint32_t>(v + 1); // bottom-right vertex
            m->indices[3 * (t) + 2] = static_cast<uint32_t>(v + 1 + V); // top-right vertex
            ++t;
            // triangle 2: v, v+1+V, v+V
            m->indices[3 * (t) + 0] = static_cast<uint32_t>(v); // bottom-left vertex
            m->indices[3 * (t) + 1] = static_cast<uint32_t>(v + 1 + V); // top-right vertex
            m->indices[3 * (t) + 2] = static_cast<uint32_t>(v + V); // top-left vertex
            ++t;
        }
    }

    // mark boundary vertices in attr (mirror MyMesh boundary array)
    // boundary order in MyMesh: 0, N, V*N, V*N+N, then for i=1..N-1: i, V*i, V*i+N, V*N+i
    
    m->attr[0] = 1.0f; // bottom-left corner
    m->attr[N] = 1.0f; // bottom-right corner
    m->attr[V * N] = 1.0f; // top-left corner
    m->attr[V * N + N] = 1.0f; // top-right corner
    m->boundary.resize(4 * N);
    size_t k = 0;
    m->boundary.data[k++] = 0;
    m->boundary.data[k++] = static_cast<uint32_t>(N);
    m->boundary.data[k++] = static_cast<uint32_t>(V * N);
    m->boundary.data[k++] = static_cast<uint32_t>(V * N + N);
    // now the rest of the boundary vertices, in the order specified above
    for (size_t i = 1; i < N; ++i) {
        m->attr[i] = 1.0f;               // bottom side
        m->attr[V * i] = 1.0f;           // right side
        m->attr[V * i + N] = 1.0f;       // left side
        m->attr[V * N + i] = 1.0f;       // top side
        m->boundary.data[k++] = static_cast<uint32_t>(i); // bottom side
        m->boundary.data[k++] = static_cast<uint32_t>(V * i); // right side
        m->boundary.data[k++] = static_cast<uint32_t>(V * i + N); // left side
        m->boundary.data[k++] = static_cast<uint32_t>(V * N + i); // top side
    }
}
