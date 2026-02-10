#pragma once

#include "vec3.h"

void inline mass_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    /* Computation of ||AB x AC|| */
    double K_360 = norm(cross(AB, AC)) / 360;

    if (K_360 == 0)
        throw std::runtime_error("The two vectors must be linearly independent.");
    else
    {
        /* Computation of the coefficients in the bloc vertex-vertex */
        M[0] = K_360 * 6;
        M[1] = K_360 * 6;
        M[2] = K_360 * 6;

        M[6] = -K_360;
        M[7] = -K_360;
        M[11] = -K_360;

        /* Computation of the coefficients in the bloc vertex-edge */
        M[8] = 0;
        M[13] = 0;
        M[17] = 0;

        M[9] = -K_360 * 4;
        M[10] = 0;
        M[14] = -K_360 * 4;

        M[12] = 0;
        M[15] = -K_360 * 4;
        M[16] = 0;

        /* Computation of the coefficients in the bloc edge-edge */
        M[3] = K_360 * 32;
        M[4] = K_360 * 32;
        M[5] = K_360 * 32;

        M[18] = K_360 * 16;
        M[19] = K_360 * 16;
        M[20] = K_360 * 16;

        /* Computation of the symmetric part */
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < i; j++)
                M[6 * i + j] = M[6 * j + i];
    }
}
