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
	/* Computation of ||AB x AC|| */
	double det = pow((AB[1] * AC[2] - AB[2] * AC[1]), 2) + pow((AB[2] * AC[0] - AB[0] * AC[2]), 2) + pow((AB[0] * AC[1] - AB[1] * AC[0]), 2);
	det = pow(det, 0.5);
	
	if (det == 0)
		throw std::runtime_error("The two vectors must be linearly independent.");
	else
	{
		/* Computation of the upper coefficients */
		S[0] = norm2(AC - AB) / (2 * det);
		S[1] = -dot(AC, AC - AB) / (2 * det);
		S[2] = dot(AB, AC - AB) / (2 * det);

		S[4] = norm2(AC) / (2 * det);
		S[5] = -dot(AB, AC) / (2 * det);

		S[8] = norm2(AB) / (2 * det);

		/* Computation of the lower coefficients */
		for (int j = 0; j < 3; j++)
			for (int i = j + 1; i < 3; i++)
				S[3 * i + j] = S[3 * j + i];
	}
}
