#pragma once

#include "vec3.h"

/* Given a triangle ABC, computes the (symmetric) 3x3 stiffness matrix S s.t.
 *
 *   S_{ij} := \int_{ABC} \nabla \phi_i \cdot \nabla \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C
 * are the shape functions of the P1 Lagrange element associated
 * to ABC.
 *
 * Input : the vectors AB and AC.
 * Output: the six coefficients S_{00} S_{11} S_{22} S_{01} S_{12} S_{20},
 *         corresponding to the interactions A<->A, B<->B, C<->C, A<->B, B<->C,
 *         C<->A
 *
 * Idea behind computation :
 * -------------------------
 *
 * Hiden for now.
 *
 */
void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
	// Squared lengths and mutual dot product
    double ABAB = norm2(AB);
    double ACAC = norm2(AC);
    double ABAC = dot(AB, AC);

    // Determinant term proportional to squared area
    double det = ABAB * ACAC - ABAC * ABAC;

    // Guard against degenerate (collinear) triangle
    if (det <= 0.0) {
        for (int i = 0; i < 6; ++i)
            S[i] = 0.0;
        return;
    }

    // Common multiplicative factor
    double mult = 0.5 / std::sqrt(det);

    ABAB *= mult;
    ACAC *= mult;
    ABAC *= mult;

    // Fill stiffness coefficients
    S[0] = ACAC - 2.0 * ABAC + ABAB;  // S00
    S[1] = ACAC;                      // S11
    S[2] = ABAB;                      // S22
    S[3] = ABAC - ACAC;               // S01
    S[4] = -ABAC;                     // S12
    S[5] = ABAC - ABAB;               // S20

}
