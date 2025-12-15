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
template <typename T>
void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    // 1. Passing from triangle in 3D to triangle in 2D
                                    // A = (0, 0)
	T x_b = norm(AB);               // B = (||AB||, 0)
    T x_c = dot(AC, AB) / x_b;      // C = (AC . AB / ||AB||, 
    //                                      ||AC - (AC . AB / ||AB||^2) * AB ||)
    T x_y = norm(AC - (x_c / norm(AB)) * AB);

    // 2. Passing from any-shape triangle in 2D to "canonic" 2D triangle:
    /*  C = (0, 1)
     *  | \
     *  |   \
     *  A ---- B = (1, 0)   and A = (0, 0)
     *
     * This transformation is performed with the following matrix:
     *          | x_b   x_c |
     * tau    = |           |
     *          |  0    y_c |
     */

    T tau_det = x_b*y_c;    // Jacobian of the trasformation is the determinant of 
                                // the matrix

    // Phi_A = -x -y -1  --> grad(Phi_A) = (-1, -1)
    // Phi_B = x         --> grad(Phi_B) = (1, 0)
    // Phi_C = y         --> grad(Phi_C) = (0, 1)

    // Finally, the integration:

    // diagonal
    M[0*3 + 0] = (1.0/12.0) * tau_det;  // int(Phi_A * Phi_A)
    M[1*3 + 1] = (1.0/12.0) * tau_det;  // int(Phi_B * Phi_B)
    M[2*3 + 2] = (1.0/12.0) * tau_det;  // int(Phi_C * Phi_C)
    // upper triangle
    M[0*3 + 1] = (1.0/24.0) * tau_det;  // int(Phi_A * Phi_B)
    M[0*3 + 2] = (1.0/24.0) * tau_det;  // int(Phi_A * Phi_C)
    M[1*3 + 2] = (1.0/24.0) * tau_det;  // int(Phi_B * Phi_C)
    // lower triangle (simmetric)
    M[1*3 + 0] = (1.0/24.0) * tau_det;  // int(Phi_B * Phi_A)
    M[2*3 + 0] = (1.0/24.0) * tau_det;  // int(Phi_C * Phi_A)
    M[2*3 + 1] = (1.0/24.0) * tau_det;  // int(Phi_B * Phi_C)
}
