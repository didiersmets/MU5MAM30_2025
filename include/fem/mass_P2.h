#pragma once

#include "vec3.h"

void inline mass_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    /* Computation of ||AB x AC|| / 360 */
    double K_360 = norm(cross(AB, AC)) / 360;

    if (K_360 == 0)
        throw std::runtime_error("The two vectors must be linearly independent.");
    else
    {
        /* Computation of the coefficients in the bloc vertex-vertex */
        M[0] = K_360 * 6;
        M[1] = -K_360;
        M[2] = -K_360;

        M[7] = K_360 * 6;
        M[8] = -K_360;

        M[14] = K_360 * 6;

        /* Computation of the coefficients in the bloc vertex-edge */
        M[3] = 0;
        M[4] = -K_360 * 4;
        M[5] = 0;

        M[9] = 0;
        M[10] = 0;
        M[11] = -K_360 * 4;

        M[15] = -K_360 * 4;
        M[16] = 0;
        M[17] = 0;

        /* Computation of the coefficients in the bloc edge-edge */
        M[21] = K_360 * 32;
        M[22] = K_360 * 16;
        M[23] = K_360 * 16;

        M[28] = K_360 * 32;
        M[29] = K_360 * 16;

        M[35] = K_360 * 32;

        /* Computation of the symmetric part */
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < i; j++)
                M[6 * i + j] = M[6 * j + i];
    }
}
