#include "mesh.h" 
#include "cube.h"  // if not already included in mesh_io.h
#include <string>

int main() {
    Mesh m;
    load_cube(m, 1);

    save_to_obj(m, "output.obj");

    return 0;
}