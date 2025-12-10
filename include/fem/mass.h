#pragma once

#include "vec3.h"

/* Solves a 2 by 2 linear system made like so:
 *     | a   b |  | x |   -  | rh1 |
 *     | c   d |  | y |   -  | rh2 |
 * using gaussian elimination
 */
template <typename T>
void inline ge_2by2(T a, T b, T c, T d, T rh1, T rh2, T* x, T* y) {
    T factor = c / a;

    T d2   = d   - factor * b;
    T rh2_ = rh2 - factor * rh1;

    // Solving for y
    T y_val = rh2_ / d2;

    // Back substitution to obtain x
    T x_val = (rh1 - b * y_val) / a;

    *x = x_val;
    *y = y_val;
}

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
	T a, b, c, d, rh1, rh2;
	


}
