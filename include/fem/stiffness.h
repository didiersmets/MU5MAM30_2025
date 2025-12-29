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

	/*
	Stiffness local matrix is defined as:
	integral over ABC of grad(phi_i) . grad(phi_j) dA

	Since basis functions are linear their gradient will just be a constant vector.
	As we computed them in the /mass.h file we can compute their gradients as:
	grad(phi_A) = grad( 1 - x - y ) = [ -1 , -1 ] and so on...
	We then compute the dot products between these gradients we integrate over the triangle
	and as before we scale by the area of the triangle divided by the area of the reference triangle (1/2)
	*/

	/* Computation of ||AB x AC|| */
	double det_2 = 2.f * norm(cross(AB, AC));
	assert(det_2 != 0);
	// diagonal 
	S[0] = norm(AC - AB) * norm(AC - AB) / det_2;
	S[1] = norm(AC) * norm(AC) / det_2;
	S[2] = norm(AB) * norm(AB) / det_2;

	// off-diagonal  ( symmetric )
	S[3] = -dot(AC - AB, AC) / det_2;
	S[4] = -dot(AB, AC) / det_2;
	S[5] = dot(AC - AB, AB) / det_2;
}
