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
	Vec3d n = cross(AB,AC);
	double aire = norm2(n);
	M[0] = (1/12)*aire; M[1] = (1/24)*aire;M[2] = (1/24)*aire;
	M[3] = (1/24)*aire;M[4] = (1/12)*aire;M[5] = (1/24)*aire;
	M[6] = (1/24)*aire;M[7] = (1/24)*aire;M[8] = (1/12)*aire;
	
}
