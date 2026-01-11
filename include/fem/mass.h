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
	// double aire = 1/2 * norm(cross(AB, AC)); // cross correspond au produit vectoriel en 3D, les fonctions norm et cross sont définies dans vec3.h
	// double coeff = aire / double(12);
	// int n = 3; 
	// for (int i = 0; i<n; i++){
	// 	for (int j = 0; j<n; j++){
	// 		if (i == j){
	// 			M[3*i + j] = 2*coeff;
	// 			}
	// 		else{
	// 			M[3*i + j] = coeff;
	// 		}
	// 	}
	// }
	double detJacobienne = norm(cross(AB,AC));

	M[0] = detJacobienne/12.0; /*Coefficients diagonaux*/
	M[1] = detJacobienne/24.0; /*Autres coefficients*/
	
}
