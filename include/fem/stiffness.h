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
	double norme2_AC = norm2(AC);
	double norme2_AB = norm2(AB);
	double pdt_scalaire = dot(AC,AB);
	double norme_BC = norme2_AC - 2 * pdt_scalaire + norme2_AB; 

	double Aire = (0.5)*(norm(cross(AB,AC)));

	/* La matrice est symetrique donc on calcule que 6 coefficients*/
    
    /* On ecrit BC = AB - AC */
    /* int (gradient phiA.gradient phiA) */
	S[0] =  norme_BC /(4*Aire);  
    /* int (gradient phiA.gradient phiB) */
	S[1] = (pdt_scalaire - norme2_AC)/(4*Aire); 
    /* int (gradient phiA.gradient phiC) */
	S[2] = (norme2_AB - pdt_scalaire )/(4*Aire); 
    /* int (gradient phiB.gradient phiB) */
	S[3] =  norme2_AC /(4*Aire) ; 
    /* int (gradient phiB.gradient phiC) */
	S[4] =	pdt_scalaire /(4*Aire); 
    /* int (gradient phiC.gradient phiC) */
	S[5] =  norme2_AB/(4*Aire);                  

}

