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
	/* Your implementation goes here ! */

	/*
	The mass matrix has a lot of redundancy. For polinomials of degree 1, over a triangle (ABC) 
	the computation is the following:

	M_ij = integral_over_triangle_ABC_of phi_i * phi_j
		 = 2 * Area(ABC) * integral_over_reference_triangle_of phi_i_hat * phi_j_hat
		 = 2 * Area(ABC) * (1/24 if i != j else 1/12)

	This means that all the diagonal values are the same, and all the off-diagonal values are the same.
	Thus, we only need to compute two values: one for the diagonal and one for the off-diagonal.
	*/

	// Diagonal entry
	const auto two_area = norm(cross(AB, AC));
	M[0] = two_area / 12.0f;

	// Off-diagonal entry
	M[1] = two_area / 24.0f;
}
