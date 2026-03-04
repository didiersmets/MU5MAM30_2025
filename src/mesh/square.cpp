#include <assert.h>
#include <stdint.h>
#include <stdio.h>



#include "math_utils.h"
#include "mesh.h"
#include "mesh/square.h"
#include "vec3.h"


static inline uint32_t vid(uint32_t i, uint32_t j, uint32_t nx) {
    // nx=number of cell in x
    return j * (nx + 1) + i;
}

int load_square_grid(Mesh &m, uint32_t nx, uint32_t ny)
{
    if (nx == 0||ny == 0) return -1;//evoiding limit cases
   
    const uint32_t nverts = (nx + 1) * (ny + 1);//number of verices
    const uint32_t ntris  = 2 * nx * ny;//number of triangles (2 per cell)

    m.positions.resize(nverts);// allocating memory for verteces 
    m.indices.resize(3u * ntris);// allocating memory for indices (3 per triangle)

    // vertex on[0, 1]x[0,1], z=0 (easier for visualization)
    for (uint32_t j = 0; j <= ny; ++j) {
         double y = (double)j / (double)ny;
    for (uint32_t i = 0; i <= nx; ++i) {
        double x = (double)i / (double)nx;
        m.positions[vid(i, j, nx)] = Vec3{ (float)x, (float)y, 0.0f };
    }
    }

    // Triangoli: ogni cella (i,j) -> due triangoli (a,b,c) e (a,c,d)
    uint32_t k=0;
    for (uint32_t j=0; j<ny; ++j) {
        for (uint32_t i = 0; i < nx; ++i) {

            //prende vertici della cella 
            uint32_t a = vid(i,   j,   nx);
            uint32_t b = vid(i+1, j,   nx);
            uint32_t c = vid(i+1, j+1, nx);
            uint32_t d = vid(i,   j+1, nx);

            // Triangolo 1
            m.indices[k++] = a;
            m.indices[k++] = b;
            m.indices[k++] = c;

            // Triangolo 2
            m.indices[k++] = a;
            m.indices[k++] = c;
            m.indices[k++] = d;
        }
    }

    return 0;
}