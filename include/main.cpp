#include <cstdio>
#include <cmath>
#include <cassert>
#include "mesh.h"

// --- Helper: Save Mesh to OBJ format ---
// This allows you to view your mesh in Blender/MeshLab
void save_obj(const Mesh& m, const char* filename) {
    FILE* f = fopen(filename, "w");
    if (!f) {
        printf("Failed to open file %s for writing!\n", filename);
        return;
    }

    // 1. Write Vertices (v x y z)
    for (size_t i = 0; i < m.vertex_count(); ++i) {
        // Access the vertex directly
        const Vec3& v = m.positions.data[i]; 
        fprintf(f, "v %f %f %f\n", v.x, v.y, v.z);
    }

    // 2. Write Faces (f i1 i2 i3)
    // Note: OBJ files are 1-indexed, so we add +1 to our indices
    for (size_t i = 0; i < m.index_count(); i += 3) {
        uint32_t i0 = m.indices.data[i + 0] + 1;
        uint32_t i1 = m.indices.data[i + 1] + 1;
        uint32_t i2 = m.indices.data[i + 2] + 1;
        fprintf(f, "f %u %u %u\n", i0, i1, i2);
    }

    fclose(f);
    printf("Saved mesh to %s\n", filename);
}

// --- Main Test ---
int main() {
    Mesh sphereMesh;
    
    // 1. Generate the Sphere
    // subdiv = 10 gives a reasonably smooth sphere (10x10 grids per face)
    printf("Generating sphere...\n");

    return 0;
}