#pragma once

#include "vec3.h"

/* Given a triangle ABC, computes the (symmetric) 3x3 mass M s.t.
 *
 *   M_{ij} := \int_{ABC} \phi_i \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C
 * are the shape functions of the P1 Lagrange element associated
 * to ABC.
 */
void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
	/* Computation of ||AB x AC|| */
	double K_12 = norm(cross(AB, AC)) / 12;

	if (K_12 == 0)
		throw std::runtime_error("The two vectors must be linearly independent.");
	else
	{
		/* Computation of the coefficients */
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				M[3 * i + j] = (i == j) ? K_12 : K_12 / 2;
	}
}
