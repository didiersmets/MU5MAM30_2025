#include "mesh.h"
#include "cube.h"
#include "sphere.h"
#include <iostream>

using namespace std;

// -----------------------------------
void load_sphere(Mesh &m, size_t subdiv)
{
    /* Load the cube */
    load_cube(m, subdiv);

    /* Normalize the cube to mesh the sphere */
    for (size_t i = 0; i < m.positions.size; i++)
        m.positions[i] /= norm(m.positions[i]);
}
