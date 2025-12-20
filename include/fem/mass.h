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
	Vec3d cross_product = cross(AB, AC);
	double triangle_area = 0.5 * norm(cross_product);

	double diag_coef = triangle_area / 6.0;
	double off_diag_coef = triangle_area / 12.0;
	M[0] = diag_coef;
	M[1] = off_diag_coef;
	M[2] = off_diag_coef;
	M[3] = off_diag_coef;
	M[4] = diag_coef;
	M[5] = off_diag_coef;
	M[6] = off_diag_coef;
	M[7] = off_diag_coef;
	M[8] = diag_coef;
}
