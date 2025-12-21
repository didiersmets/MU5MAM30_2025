#include <iostream>
// #include "include/mesh/cube.h"
// #include "include/mesh/sphere.h"
// #include "include/mesh/mesh.h"
// #include "include/fem/P1.h"

#include "include/matrix/sparse_matrix.h"
#include "include/fem/P1.h"
#include "include/mesh/cube.h"
#include "include/mesh/sphere.h"
using namespace std;

int main()
{
    size_t subdiv = 2;

    Mesh cube;
    load_cube(cube, subdiv);

    bool print_cube = true;
    bool print_sphere = false;

    // cout << cube.vertex_count() << endl;

    if (print_cube)
    {
        cout << "Cube positions : " << endl;
        for (size_t i = 0; i < cube.positions.size; i++)
            cout << cube.positions[i].x << " " << cube.positions[i].y << " " << cube.positions[i].z << endl;

        cout << endl;

        cout << "Cube indices : " << endl;
        for (size_t i = 0; i < floor(cube.indices.size / 3); i++)
            cout << cube.indices[3 * i] << " " << cube.indices[3 * i + 1] << " " << cube.indices[3 * i + 2] << endl;
    }

    if (print_sphere)
    {
        Mesh sphere;
        load_sphere(sphere, subdiv);

        cout << endl;
        cout << "Sphere : " << endl;
        for (size_t i = 0; i < sphere.positions.size; i++)
            cout << sphere.positions[i].x << " " << sphere.positions[i].y << " " << sphere.positions[i].z << endl;
    }
    return 0;
}