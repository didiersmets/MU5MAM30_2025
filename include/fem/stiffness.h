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
	/* Your implementation goes here */
    const double ab2 = dot(AB, AB);
    const double ac2 = dot(AC, AC);
    const double ab_ac = dot(AB, AC);

    const double jac = std::sqrt(ab2 * ac2 - ab_ac * ab_ac);
    const double inv2A = 1.0 / jac;   
    const double coeff = 0.5 * inv2A; 

    S[0] = coeff * (ab2 - 2.0 * ab_ac + ac2);
    S[1] = coeff * ac2;
    S[2] = coeff * ab2;
    S[3] = coeff * (ab_ac - ac2);
    S[4] = coeff * (-ab_ac);
    S[5] = coeff * (ab_ac - ab2);           
}
