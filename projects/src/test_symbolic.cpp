#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "mesh.h"
#include "mesh_bounds.h"
#include "mesh_gpu.h"
#include "mesh_io.h"
#include "cube.h"
#include "symbolic.h"
#include "sparse_matrix.h"
#include "array.h"
#include "P1.h"

#include <unordered_set>

//print sparisty pattern by ChatGPT
void printSparsityPattern(const CSRPattern& A, char filled='*', char empty=' ')
{
    for (size_t i = 0; i < A.rows; ++i)
    {
        uint32_t k = A.row_start[i];
        uint32_t k_end = A.row_start[i + 1];

        for (size_t j = 0; j < A.cols; ++j)
        {
            if (k < k_end && A.col[k] == j)
            {
                std::cout << filled;
                ++k;
            }
            else
            {
                std::cout << empty;
            }
        }
        std::cout << '\n';
    }
}


int main(int argc, char** argv) {
    bool print_mode = false;
    size_t subdiv = 2;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-p") == 0 ||
            std::strcmp(argv[i], "--print") == 0) {
            print_mode = true;
        }
        else if (std::strcmp(argv[i], "-s") == 0 ||
                 std::strcmp(argv[i], "--subdiv") == 0) {
            if (i + 1 >= argc) {
                std::cerr << "Error: --subdiv requires a value\n";
                return 1;
            }
            subdiv = static_cast<size_t>(std::atoi(argv[++i]));
        }
        else {
            std::cerr << "Unknown argument: " << argv[i] << "\n";
            return 1;
        }
    }

    Mesh mesh, mesh_disect;
    load_cube(mesh, subdiv);
    load_cube_nested_dissect(mesh_disect, subdiv);

    CSRPattern P, P_D,  L, L_D;
    build_P1_CSRPattern(mesh, P);
    build_P1_CSRPattern(mesh_disect, P_D);

    TArray<uint32_t> parent(P.rows);
    construct_etree(P, parent);
    construct_L_sparsity_pattern(P, L, parent);

    TArray<uint32_t> parent_D(P_D.rows);
    construct_etree(P_D, parent_D);
    construct_L_sparsity_pattern(P_D, L_D, parent_D);

    if (print_mode && P.rows > 50) {
    std::cerr << "Matrix too large to print safely\n";
    return 1;
    }

    if (print_mode) {
        std::cout << "Input matrix P:\n";
        printSparsityPattern(P);

        std::cout << "Elimination tree (parent):\n";
        for (size_t i = 0; i < parent.size; ++i) {
            std::cout << parent[i] << " ";
        }
        std::cout << "\n";

        std::cout << "Cholesky sparsity pattern L:\n";
        printSparsityPattern(L);
    } else {
        std::cout << "subdiv:  " << subdiv << "\n";
        std::cout << "n: " << P.rows << "\n";
        std::cout << "P.nnz = " << P.nnz << "\n";
        std::cout << "L.nnz = " << L.nnz << "\n";
        std::cout << "L_D.nnz = " << L_D.nnz << "\n";
    }

    return 0;
}


