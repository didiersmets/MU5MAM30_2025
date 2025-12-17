#include "cube.h"
#include "duplicate_verts.h"
#include <cmath>

static void add_face(Mesh &m, size_t subdiv, Vec3 origin, Vec3 u_dir, Vec3 v_dir) {
    size_t start_idx = m.vertex_count();

    // 1. Generate Vertices
    for (size_t y = 0; y <= subdiv; ++y) {
        for (size_t x = 0; x <= subdiv; ++x) {
            float u = (float)x / subdiv;
            float v = (float)y / subdiv;

            Vec3 p;
            p.x = origin.x + (u_dir.x * u * 2.0f) + (v_dir.x * v * 2.0f);
            p.y = origin.y + (u_dir.y * u * 2.0f) + (v_dir.y * v * 2.0f);
            p.z = origin.z + (u_dir.z * u * 2.0f) + (v_dir.z * v * 2.0f);

            m.positions.push_back(p); 
        }
    }

    // 2. Generate Indices
    size_t row_len = subdiv + 1;
    for (size_t y = 0; y < subdiv; ++y) {
        for (size_t x = 0; x < subdiv; ++x) {
            uint32_t bl = start_idx + (y * row_len) + x;
            uint32_t br = start_idx + (y * row_len) + x + 1;
            uint32_t tl = start_idx + ((y + 1) * row_len) + x;
            uint32_t tr = start_idx + ((y + 1) * row_len) + x + 1;

            m.indices.push_back(bl);
            m.indices.push_back(br);            
            m.indices.push_back(tr);

            m.indices.push_back(bl);
            m.indices.push_back(tr);
            m.indices.push_back(tl);
        }
    }
}

int load_cube(Mesh &m, size_t subdiv) {
    // If your TArray has a .clear() function, call it here:
    // m.positions.clear();
    // m.indices.clear();

    add_face(m, subdiv, {-1, -1,  1}, { 1, 0, 0}, { 0, 1, 0});
    add_face(m, subdiv, { 1, -1, -1}, {-1, 0, 0}, { 0, 1, 0});
    add_face(m, subdiv, { 1, -1,  1}, { 0, 0,-1}, { 0, 1, 0});
    add_face(m, subdiv, {-1, -1, -1}, { 0, 0, 1}, { 0, 1, 0});
    add_face(m, subdiv, {-1,  1,  1}, { 1, 0, 0}, { 0, 0,-1});
    add_face(m, subdiv, {-1, -1, -1}, { 1, 0, 0}, { 0, 0, 1});

    remove_duplicate_vertices(m);

    return 0;
}