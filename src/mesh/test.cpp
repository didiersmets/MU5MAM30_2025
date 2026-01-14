#include <iostream>
#include "mesh.h"
#include "adjacency.h"

int main() {
    Mesh m;

    // Creamos una malla simple con 4 vértices y 2 triángulos
    m.positions.resize(4);
    m.positions[0] = {0,0,0};
    m.positions[1] = {1,0,0};
    m.positions[2] = {0,1,0};
    m.positions[3] = {1,1,0};

    // Triángulos: (0,1,2) y (1,3,2)
    m.indices.resize(6);
    m.indices[0] = 0; m.indices[1] = 1; m.indices[2] = 2;
    m.indices[3] = 1; m.indices[4] = 3; m.indices[5] = 2;

    VTAdjacency adj(m);

    // Imprimir degree
    std::cout << "degree:\n";
    for (size_t i = 0; i < m.vertex_count(); i++)
        std::cout << "v" << i << ": " << adj.degree[i] << "\n";

    // Imprimir offset
    std::cout << "\noffset:\n";
    for (size_t i = 0; i < m.vertex_count(); i++)
        std::cout << "v" << i << ": " << adj.offset[i] << "\n";

    // Imprimir vtri
    std::cout << "\nvtri:\n";
    for (size_t i = 0; i < adj.vtri.size; i++)
        std::cout << i << ": (" << adj.vtri[i].next << ", " << adj.vtri[i].prev << ")\n";

    return 0;
}
