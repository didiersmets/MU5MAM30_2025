#pragma once

#include "vec3.h"
#include "sys_utils.h"

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
	ASSERT(S!= nullptr);

	Vec3d edge_opposite_to_A = AC - AB; // BC

    double area2 = 0.5 * norm(cross(AB, AC));
	constexpr double epsilon = 1e-14;
	ASSERT_ALWAYS(area2 > epsilon);

    double coefficient = 1.0 / (4.0 * area2);

	S[0] = - coefficient * dot(AC, edge_opposite_to_A); 
	S[1] = - coefficient * dot(AB, edge_opposite_to_A);                                   // S_BB
	S[2] = - coefficient * dot(AB, AC);

}
