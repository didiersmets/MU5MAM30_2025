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
    double ab2   = norm2(AB);
    double ac2   = norm2(AC);
    double ab_ac = dot(AB, AC);

    double det = ab2 * ac2 - ab_ac * ab_ac;
    double scale = 0.5 / sqrt(det);

    ab2   *= scale;
    ac2   *= scale;
    ab_ac *= scale;

    S[0] = ac2 - 2.0 * ab_ac + ab2;
    S[1] = ac2;
    S[2] = ab2;
    S[3] = ab_ac - ac2;
    S[4] = -ab_ac;
    S[5] = ab_ac - ab2;
}
