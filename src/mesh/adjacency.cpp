#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
    const uint32_t n_vertices  = static_cast<uint32_t>(m.vertex_count()); // cuántos vértices hay
    const uint32_t n_triangles = static_cast<uint32_t>(m.triangle_count()); // cuántos grupos de 3 índices hay

    // 1) degree: cuántos triángulos tocan a cada vértice

	//crea el array degree y lo llena con 0s
    degree.resize(n_vertices);
    for (uint32_t i = 0; i < n_vertices; ++i) {
        degree[i] = 0;
    }

	 //recorre todos los triángulos y cuenta el número de incidencias
    for (uint32_t t = 0; t < n_triangles; ++t) {
        uint32_t a = m.indices[3*t + 0];
        uint32_t b = m.indices[3*t + 1];
        uint32_t c = m.indices[3*t + 2];

        degree[a] += 1;
        degree[b] += 1;
        degree[c] += 1;
    }

    // 2) offset: prefijo acumulado de degree
    offset.resize(n_vertices);
    uint32_t acc = 0;
    for (uint32_t i = 0; i < n_vertices; ++i) {
        offset[i] = acc;
        acc += degree[i];
    }

    // 3) vtri: tamaño total = suma de degree = 3 * n_triangles
    vtri.resize(acc);

    // 4) cursor: para ir llenando vtri dentro del bloque de cada vértice
    TArray<uint32_t> cursor;
    cursor.resize(n_vertices);
    for (uint32_t i = 0; i < n_vertices; ++i) {
        cursor[i] = offset[i];
    }

    // 5) rellenar vtri
    for (uint32_t t = 0; t < n_triangles; ++t) {
        uint32_t a = m.indices[3*t + 0];
        uint32_t b = m.indices[3*t + 1];
        uint32_t c = m.indices[3*t + 2];

        // para el vértice a
        {
            uint32_t j = cursor[a]++;
            vtri[j].next = b;
            vtri[j].prev = c;
        }
        // para el vértice b
        {
            uint32_t j = cursor[b]++;
            vtri[j].next = c;
            vtri[j].prev = a;
        }
        // para el vértice c
        {
            uint32_t j = cursor[c]++;
            vtri[j].next = a;
            vtri[j].prev = b;
        }
    }
}
