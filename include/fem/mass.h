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

	M[0]= triangle_area / 6.0;
	M[1] = triangle_area / 12.0;
}

void inline mass_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
	Vec3d cross_product = cross(AB, AC);
	double triangle_area = 0.5 * norm(cross_product);

	M[0] = triangle_area * (6.0/180.0);
	M[1] = triangle_area * (-1.0/180.0);
	M[2] = triangle_area * (32.0/180.0);
	M[3] = triangle_area * (16.0/180.0);
	M[4] = triangle_area * (-4.0/180.0);
}
