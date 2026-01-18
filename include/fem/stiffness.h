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
	Vec3d *T {{0.,0.,0.}, AB, AC};
	float area {norm(cross(AB, AC)) /2};


	S[0] = (T[2]-T[1])*(T[2]-T[1])/(4*area);  	// S_{00}
	S[1] = (T[0]-T[2])*(T[0]-T[2])/(4*area);	// S_{11}
	S[2] = (T[1]-T[0])*(T[1]-T[0])/(4*area);	// S_{22}
	S[3] = (T[2]-T[1])*(T[0]-T[2])/(4*area);	// S_{01}
	S[4] = (T[0]-T[2])*(T[1]-T[0])/(4*area);	// S_{12}
	S[5] = (T[1]-T[0])*(T[2]-T[1])/(4*area);	// S_{20}
}
