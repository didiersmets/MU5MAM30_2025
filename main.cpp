#include "utils_debug.h"
#include "sphere.h"
#include "mesh.h"

int main() {
    Mesh mesh;
    size_t subdiv = 4;   
    size_t row_size = subdiv + 1; 

    load_sphere(mesh, subdiv);

    // print_mesh_info(mesh);
    print_square_mesh(mesh, 4, row_size);

}