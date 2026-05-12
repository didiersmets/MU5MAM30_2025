#include <iostream>
#include <string.h>
#include <functional>

#include "minimal_graph.h"
#include "my_mesh.h"
#include "tiny_blas.h"
#include "P1.h"

////////////////////////////////////////////////
//   Test function to be changed as needed    //
////////////////////////////////////////////////

static double test_f(const Vec2d &pos)
{
    double x = pos.x;
    double y = pos.y;
    return 1/3.0 * (pow(x  - 2,3) +3*(x - 2)*pow(y - 2,2));
}


int main(int argc, char **argv){

    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <mesh type> <N> <size1> <size2>" << std::endl;
        std::cerr << "Mesh type: disk or square" << std::endl;

        return 1;
    }

    const char* mesh_type = argv[1];
    int N = atoi(argv[2]);
    double size1 = 1.0, size2 = 1.0;
    if (argc > 3) {
        size1 = atof(argv[3]);
        if (size1 <= 0) {
            std::cerr << "Size1 must be a positive number." << std::endl;
            return 1;
        }
    }

    if (argc == 5) {
        size2 = atof(argv[4]);
        if (size2 <= 0) {
            std::cerr << "Size2 must be a positive number." << std::endl;
            return 1;
        }
    }

    if (N <= 0) {
        std::cerr << "Number of vertices must be a positive integer." << std::endl;
        return 1;
    }

    Mesh m;
    if (strcmp(mesh_type, "disk") == 0) {
        build_disk_mesh(&m, N, size1);
    } else if (strcmp(mesh_type, "square") == 0) {
        build_square_mesh(&m, N, size1, size2);
    } else {
        std::cerr << "Invalid mesh type. Use 'disk' or 'square'." << std::endl;
        return 1;
    }

    printf("Number of DOF : %ld\n", m.vtx_count);
    
    MinimalGraphSolver solver(m, test_f);

    solver.do_iterate_Newton(200,1e-12,0.1);
    solver.do_iterate_Picardi();    

    
    return 0;

}