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
	double Area = 0.5 * norm2(cross(AB,AC));
	double C1= Area / 6;
	double C2= Area /12;

	M[0]=C1;
	M[1]=C2;
	M[2]=C2;
	M[3]=C2;
	M[4]=C1;
	M[5]=C2;
	M[6]=C2;
	M[7]=C2;
	M[8]=C1;


}
