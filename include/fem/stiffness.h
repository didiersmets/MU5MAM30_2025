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
	double denom = 2.0f * norm(cross(AB, AC));
	const Vec3d BC = AC - AB;

	// following the notes, the stiffness matrix coefficients are given by integrals over the triangle
	// of the gradients of the shape functions. For 1st degree polynomials these are constants.
	// The computation is the following:

	// diagonal terms
	S[0] = norm(BC) * norm(BC) / denom;
	S[1] = norm(AC) * norm(AC) / denom;
	S[2] = norm(AB) * norm(AB) / denom;

	// off-diagonal terms
	S[3] = -dot(BC, AC) / denom;
	S[4] = -dot(AB, AC) / denom;
	S[5] = dot(BC, AB) / denom;
}
