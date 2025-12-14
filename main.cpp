#include <iostream>
#include "utils_debug.h"
#include "sphere.h"
#include "cube.h"
#include "mesh.h"

#include <chrono>

int main() {
    Mesh mesh;
    size_t subdiv = 16;   

    auto start = std::chrono::system_clock::now();
    load_sphere(mesh, subdiv);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time taken: " << elapsed.count() << "ms" << std::endl;

    // print_mesh_info(mesh);
    save_mesh_stl(mesh, "output_sphere.stl");

}