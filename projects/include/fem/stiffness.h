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
	Vec3d n = cross(AB,AC);
	double aire = norm2(n);
	Vec3d BC = AC - AB;
	S[0] = dot(BC,BC)/4*aire;S[1] = dot(AC,AC)/4*aire;S[2] = dot(AB,AB)/4*aire;
	S[3] = dot(-AC,BC)/4*aire;S[4] = dot(-AC,AB)/4*aire;S[5] = dot(AB,BC)/4*aire;

}
