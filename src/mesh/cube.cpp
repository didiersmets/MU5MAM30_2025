#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

int load_cube(Mesh &m, size_t subdiv)
{
    /* Check subdiv is reasonable */
    if (subdiv == 0 || subdiv > (1 << 14)) {
        return -1;
    }

    const size_t n = subdiv + 1;
    const float step = 2.0f / subdiv;

    /* Reserve memory */
    m.positions.resize(6 * n * n);
    m.indices.resize(6 * subdiv * subdiv * 6);

    size_t vbase = 0;  /* base vertex index for each face */
    size_t ibase = 0;  /* base index offset */

    auto emit_face = [&](Vec3 origin, Vec3 du, Vec3 dv)
    {
        /* Generate vertices */
        for (size_t j = 0; j < n; ++j) {
            for (size_t i = 0; i < n; ++i) {
                m.positions[vbase + j * n + i] =
                    origin
                    + du * (float(i) * step)
                    + dv * (float(j) * step);
            }
        }

        /* Generate indices */
        for (size_t j = 0; j < subdiv; ++j) {
            for (size_t i = 0; i < subdiv; ++i) {
                uint32_t v0 = uint32_t(vbase + j * n + i);
                uint32_t v1 = uint32_t(vbase + j * n + i + 1);
                uint32_t v2 = uint32_t(vbase + (j + 1) * n + i);
                uint32_t v3 = uint32_t(vbase + (j + 1) * n + i + 1);

                /* Triangle 1 */
                m.indices[ibase++] = v0;
                m.indices[ibase++] = v1;
                m.indices[ibase++] = v2;

                /* Triangle 2 */
                m.indices[ibase++] = v1;
                m.indices[ibase++] = v3;
                m.indices[ibase++] = v2;
            }
        }

        vbase += n * n;
    };

    /* +X face */
    emit_face(
        Vec3{ 1, -1, -1 },
        Vec3{ 0,  0,  1 },
        Vec3{ 0,  1,  0 }
    );

    /* -X face */
    emit_face(
        Vec3{ -1, -1,  1 },
        Vec3{ 0,  0, -1 },
        Vec3{ 0,  1,  0 }
    );

    /* +Y face */
    emit_face(
        Vec3{ -1, 1, -1 },
        Vec3{ 1,  0,  0 },
        Vec3{ 0,  0,  1 }
    );

    /* -Y face */
    emit_face(
        Vec3{ -1, -1,  1 },
        Vec3{ 1,  0,  0 },
        Vec3{ 0,  0, -1 }
    );

    /* +Z face */
    emit_face(
        Vec3{ -1, -1, 1 },
        Vec3{ 1,  0, 0 },
        Vec3{ 0,  1, 0 }
    );

    /* -Z face */
    emit_face(
        Vec3{ 1, -1, -1 },
        Vec3{ -1, 0,  0 },
        Vec3{ 0,  1,  0 }
    );

    /* Merge coincident vertices */
    remove_duplicate_vertices(m);

    return 0;
}