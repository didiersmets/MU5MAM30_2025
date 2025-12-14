#include "sphere.h"
#include "cube.h"

int load_sphere(Mesh &m, size_t subdiv) {
    load_overlapping_cube<float, true>(m, subdiv);
    remove_duplicate_vertices(m);
    return 1;
}