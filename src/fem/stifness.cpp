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
	Vec3d BC = AB - AC;
	double area = norm(cross(AB, AC))/2;

	//First store diagonal coefficients
	S[0] = norm(BC)/(4*area);
	S[1] = norm(AC)/(4*area);
	S[2] = norm(AB)/(4*area);
	
	// Then store other coefficients in order AB, AC, BC
	S[3] = dot(-1*AC, AC - AB )/(4*area);
	S[4] = dot(B, AC - AB)/ (4*area);
	S[5] = dot(-1*AC, AB)/ (4*area);
}	