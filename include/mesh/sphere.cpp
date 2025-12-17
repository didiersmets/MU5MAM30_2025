#include "sphere.h"
#include "cube.h"
#include <cmath>

int load_sphere(Mesh &m, size_t subdiv) {
    // 1. Build a cube first
    // This gives us the correct "Grid" topology wrapped around a center point
    load_cube(m, subdiv);

    // 2. Inflate the cube
    // Iterate over every vertex and normalize it (make length = 1)
    for (size_t i = 0; i < m.positions.size; ++i) {
        Vec3& p = m.positions[i];

        // Calculate current distance from center
        float len = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);

        // Project onto unit sphere
        if (len > 1e-6f) {
            p.x /= len;
            p.y /= len;
            p.z /= len;
        }
    }

    return 0;
}