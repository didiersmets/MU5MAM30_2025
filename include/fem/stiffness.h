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
 */
void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
	/* Computation of ||AB x AC|| */
	double K2 = 2 * norm(cross(AB, AC));

	if (K2 == 0)
		throw std::runtime_error("The two vectors must be linearly independent.");
	else
	{
		Vec3d BC = AC - AB;
		S[0] = norm2(BC) / (K2);	// A <-> A
		S[1] = norm2(AC) / (K2);	// B <-> B
		S[2] = norm2(AB) / (K2);	// C <-> C
		S[3] = -dot(AC, BC) / (K2); // A <-> B
		S[4] = -dot(AB, AC) / (K2); // B <-> C
		S[5] = dot(AB, BC) / (K2);	// C <-> A
	}
}
