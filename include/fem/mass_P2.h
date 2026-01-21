#pragma once

#include "vec3.h"

/* Given a triangle ABC, whose middle edges points are counterclockwise given by D, E and F,
 * computes the (symmetric) 6x6 mass M s.t.
 *
 *   M_{ij} := \int_{ABC} \phi_i \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C,
 * \phi_3 := \phi_D, \phi_3 := \phi_E and \phi_5 := \phi_F
 * are the shape functions of the P2 Lagrange element associated
 * to ABC.
 */
void inline mass_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    /* Computation of ||AB x AC|| */
    double K = norm(cross(AB, AC));

    if (K == 0)
        throw std::runtime_error("The two vectors must be linearly independent.");
    else
    {
        /* Computation of the coefficients in the bloc (0,0) */
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                M[3 * i + j] = (i == j) ? 6 * K / 360 : -K / 360;

        /* Computation of the coefficients in the bloc (0,1) */
        for (int i = 0; i < 3; i++)
            for (int j = 3; j < 6; j++)
                M[3 * i + j] = (i == ((j - 4) % 3 + 3) % 3) ? -4 * K / 360 : 0;

        /* Computation of the coefficients in the bloc (1,0) */
        for (int i = 3; i < 6; i++)
            for (int j = 0; j < 3; j++)
                M[3 * i + j] = M[3 * j + i];

        /* Computation of the coefficients in the bloc (1,1) */
        for (int i = 3; i < 6; i++)
            for (int j = 3; j < 6; j++)
                M[3 * i + j] = (i == j) ? 32 * K / 360 : 16 * K / 360;
    }
}
