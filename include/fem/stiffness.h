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
	double ABAB = norm2(AB);
	double ACAC = norm2(AC);
	double ABAC = dot(AB, AC);
	double mult = 0.5 / sqrt(ABAB * ACAC - ABAC * ABAC);
	ABAB *= mult;
	ACAC *= mult;
	ABAC *= mult;

	S[0] = ACAC - 2 * ABAC + ABAB;
	S[1] = ACAC;
	S[2] = ABAB;
	S[3] = ABAC - ACAC;
	/* Note the chosen order : (B,C)-> 4 and (C,A) -> 5 */
	S[4] = -ABAC;
	S[5] = ABAC - ABAB;
}

void inline stiffness_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
	Vec3d cross_product = cross(AB, AC);
	double cp_area = norm(cross_product);
	double trig_area = cp_area * 0.5;
	Vec3d n = cross_product * (1/cp_area);

	Vec3d BC = AC - AB;
	Vec3d CA = -1.0 * AC;

	Vec3d grad_L1 = cross(n, BC) * (1/cp_area);
	Vec3d grad_L2 = cross(n, CA) * (1/cp_area);
	Vec3d grad_L3 = cross(n, AB) * (1/cp_area);

	double dot_12 = dot(grad_L1, grad_L2);
	double dot_13 = dot(grad_L1, grad_L3);
	double dot_23 = dot(grad_L2, grad_L3);

	double dot_11 = norm2(grad_L1);
	double dot_22 = norm2(grad_L2);
	double dot_33 = norm2(grad_L3);

	// Diag vertex
	S[0] = dot_11 * trig_area;
	S[1] = dot_22 * trig_area;
	S[2] = dot_33 * trig_area;

	// Diag mid-point
	S[3] = 8 * (dot_11 + dot_12 + dot_22) * trig_area / 3.0;
	S[4] = 8 * (dot_22 + dot_23 + dot_33) * trig_area / 3.0;
	S[5] = 8 * (dot_11 + dot_13 + dot_33) * trig_area / 3.0;

	// Upper trig vertex
	S[6] = -1.0 * dot_12 * trig_area / 3.0;
	S[7] = -1.0 * dot_23 * trig_area / 3.0;
	S[8] = -1.0 * dot_13 * trig_area / 3.0;

	// Upper trig mid-point
	S[9] = 8 * dot_13 * trig_area / 3.0;
	S[10] = 8 * dot_12 * trig_area / 3.0;
	S[11] = 8 * dot_23 * trig_area / 3.0;

	// Upper block
	S[12] = 4 * dot_12 * trig_area / 3.0;
	S[13] = 4 * dot_23 * trig_area / 3.0;
	S[14] = 4 * dot_13 * trig_area / 3.0;
}