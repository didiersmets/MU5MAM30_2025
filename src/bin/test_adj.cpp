#include "adjacency.h"
#include "mesh.h"
#include <iostream>

int main() {
    Mesh m;
    //                   A              B              C
    m.positions.push_back({0.5,0.5,2.0});
    m.positions.push_back({0.0,0.0,1.0});
    m.positions.push_back({1.0,0.0,1.0});
    m.positions.push_back({1.0,1.0,1.0});
    m.positions.push_back({0.0,1.0,1.0});
    m.positions.push_back({0.0,0.0,0.0});
    m.positions.push_back({1.0,0.0,0.0});
    m.positions.push_back({1.0,1.0,0.0});
    m.positions.push_back({0.0,1.0,0.0});
    int data[] = {
        0,1,2, 0,2,3, 0,3,4, 0,4,1,
        4,5,1, 4,8,5, 1,5,6, 1,6,2,
        2,6,7, 2,7,3, 3,7,8, 3,8,4,
        5,7,6, 5,8,7
    };

    for (int v : data)
        m.indices.push_back(v);

    // TArray<uint32_t> exp_degrees = {4,5,5,5,5,5,4,5,4};
    // TArray<uint32_t> exp_offsets = {0,4,9,14,19,24,29,33,38};
	// struct VTri {
	// 	uint32_t next;
	// 	uint32_t prev;
	// };
    // Tarray<uint32_t> exp_vtri = { {1,2}, {2,3}, {3,4}, {4,1},
    //                               {2,0}, {0,4}, {4,5}, {5,6}, {6,2}
    //                               {0,1}, {3,0}, {1,6}, {6,7}, {7,3},
    //                               {0,2}, {4,0}, {2,7}, {7,8}, {8,4},
    //                               {0,3}, {1,0}, {8,5}, {5,1}, {3,8},
    //                               {4,8}, {1,4}, {6,1}, {6,7}, {8,7},
    //                               {1,5}, {2,1}, {7,2}, {5,7},
    //                               {2,6}, {3,2}, {7,3}, {6,5}, {5,8},
    //                               {5,4}, {3,7}, {4,3}, {7,5}
    //                              };

    VTAdjacency adj(m);

    std::cout << "Degrees: " << std::endl;
    std::cout << "-- Expected: 4 5 5 5 5 5 4 5 4" << std::endl;
    std::cout << "-- Got:      ";
    for (int i=0; i<9; i++) {
        std::cout << adj.degree[i] << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "Offsets: " << std::endl;
    std::cout << "-- Expected: 0 4 9 14 19 24 29 33 38" << std::endl;
    std::cout << "-- Got:      ";
    for (int i=0; i<9; i++) {
        std::cout << adj.offset[i] << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "Vtri size: " << std::endl; 
    std::cout << "-- Expected: 42" << std::endl; 
    std::cout << "-- Got:        " << adj.vtri.size << std::endl << std::endl;
    std::cout << "Vtri: ";
    for (int i=0; i<adj.vtri.size; i++) {
        std::cout << adj.vtri[i].next << "," << adj.vtri[i].prev << " ";
    }
    std::cout << std::endl;

}