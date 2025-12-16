#pragma once

#include "vec3.h"

/* Given a triangle ABC, computes the (symmetric) 3x3 mass M s.t.
 *
 *   M_{ij} := \int_{ABC} \phi_i \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C
 * are the shape functions of the P1 Lagrange element associated
 * to ABC.
 *
 * Idea behind computation :
 * -------------------------
 *
 * Hiden for now.
 */


void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    // 1) Área del triángulo ABC
    // Área = 0.5 * || AB x AC ||
    const double area = 0.5 * norm(cross(AB, AC));

    // 2) Coeficientes de la matriz de masa P1
    const double diag = area / 6.0;   // M_ii
    const double off  = area / 12.0;  // M_ij, i != j

    // 3) Relleno de la matriz (row-major)
    M[0] = diag;  M[1] = off;   M[2] = off;
    M[3] = off;   M[4] = diag;  M[5] = off;
    M[6] = off;   M[7] = off;   M[8] = diag;
}
