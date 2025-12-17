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
    double area2 = 2*norm(cross(AB,AC));

    double AB_2 = norm2(AB);
    double AC_2 = norm2(AC);
    double BC_2 = AC_2 + AB_2 - 2*dot(AB,AC);
    double AB_AC = dot(AB,AC);

    S[0] = BC_2 / area2;
    S[1] = AC_2 / area2;
    S[2] = AB_2 / area2;
    S[3] = -AB_AC / area2;
    S[4] = -AC_2 + AB_AC / area2;
    S[5] = -AB_2 + AB_AC / area2;
}
