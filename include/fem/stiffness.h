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
	/* Your implementation goes here */
	double dot_ABAC = dot(AB,AC);
	double dot_ABAB = norm2(AB);
	double dot_ACAC = norm2(AC);

	Vec3d cross_product = cross(AB,AC);
	double cp_area = norm(cross_product);

	double scaling = -0.5 / cp_area;

	double s_BC = scaling * dot_ABAC;
	double s_CA = scaling * (dot_ABAB - dot_ABAC);
	double s_AB = scaling * (dot_ACAC - dot_ABAC);

	S[0] = -(s_AB + s_CA);
	S[1] = -(s_AB + s_BC);
	S[2] = -(s_BC + s_CA);

	S[3] = s_AB;
	S[4] = s_BC;
	S[5] = s_CA;
}
