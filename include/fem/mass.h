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
	double area = norm(cross(AB, AC)) * 0.5;
	double diag = area / 12.0;
	double non_diag = area / 24.0;

	M[0] = diag;
	M[1] = non_diag;
}