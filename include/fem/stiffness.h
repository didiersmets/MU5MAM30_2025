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
	// cross = AB x AC
	Vec3d BC = AC - AB;

	Vec3d cross;
	cross.x = AB.y * AC.z - AB.z * AC.y;
	cross.y = AB.z * AC.x - AB.x * AC.z;
	cross.z = AB.x * AC.y - AB.y * AC.x;

	double cross_norm = std::sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
	double area = 0.5 * cross_norm; //l'aire du triangle ABC vaut 1/2 |AB x AC|

	double bcbc = BC.x * BC.x + BC.y * BC.y + BC.z * BC.z;	
	double acac = AC.x * AC.x + AC.y * AC.y + AC.z * AC.z;
	double abab = AB.x * AB.x + AB.y * AB.y + AB.z * AB.z;
	double abac = AB.x * AC.x + AB.y * AC.y + AB.z * AC.z;
	double abbc = AB.x * BC.x + AB.y * BC.y + AB.z * BC.z;
	double acbc = AC.x * BC.x + AC.y * BC.y + AC.z * BC.z;

	double coeff = 1.0/(4.0*area) ;

	S[0] = coeff * bcbc;  S[1] = coeff * acac;  S[2] = coeff * abab;
	S[3] = -coeff * acbc; S[4] = -coeff * abac; S[5] = -coeff * abbc;
    
}



