#include "cube.h"



int load_cube(Mesh &m, size_t subdiv) {
    load_overlapping_cube<float, false>(m, subdiv);
    remove_duplicate_vertices(m);
    return 1;
}
