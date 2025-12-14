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
     *       | 1/x_b   -x_c/(x_b*y_c) |                                  | x_b   x_c |
     * tau = |                        |   obtained by inverting tau^-1 = |           |
     *       |  0         1/y_c       |                                  |  0    y_c |
     */

    T tau_det = 1/(x_b*y_c);    // Jacobian of the trasformation is the determinant of 
                                // the matrix

    // Phi_A = -x -y -1  --> grad(Phi_A) = (-1, -1)
    // Phi_B = x         --> grad(Phi_B) = (1, 0)
    // Phi_C = y         --> grad(Phi_C) = (0, 1)

    // Finally, the integration:              | (x_b + x_c)^2  x_c*y_c |   | a     b |
	// We need to build tau^(-1) * tau^(-T) = |                        | = |         |
	//										  |    x_c*y_c       y_c^2 |   | b     c |
	T a = (x_b + x_c) * (x_b + x_c);
	T b = x_c * y_c;
	T c = y_c * y_c;

    // diagonal
    S[0*3 + 0] = (a+b+c) * tau_det / 2;  // int(grad(Phi_A) * grad(Phi_A))
    S[1*3 + 1] = a       * tau_det / 2;  // int(grad(Phi_B) * grad(Phi_B))
    S[2*3 + 2] = c       * tau_det / 2;  // int(grad(Phi_C) * grad(Phi_C))
    // upper triangle
    S[0*3 + 1] = (-1) * (a+b) * tau_det / 2;  // int(grad(Phi_A) * grad(Phi_B))
    S[0*3 + 2] = (-1) * (b+c) * tau_det / 2;  // int(grad(Phi_A) * grad(Phi_C))
    S[1*3 + 2] = b            * tau_det / 2;  // int(grad(Phi_B) * grad(Phi_C))
    // lower triangle (simmetric)
    S[1*3 + 0] = (-1) * (a+b) * tau_det / 2;  // int(grad(Phi_B) * grad(Phi_A))
    S[2*3 + 0] = (-1) * (b+c) * tau_det / 2;  // int(grad(Phi_C) * grad(Phi_A))
    S[2*3 + 1] = b            * tau_det / 2;  // int(grad(Phi_B) * grad(Phi_C))
}
