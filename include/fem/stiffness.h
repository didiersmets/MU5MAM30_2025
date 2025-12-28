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
	// attention, on appelle une matrice S qui est une FEMatrix et qui est de type P1_sym
	// pourquoi : on a que la matrice stiffness est une matrice symétrique car $int \nabla \phi_A \phi_B = int \nmbla \phiB \phi_A$
	//
	//   A	B	C
	//A		0	1
	//B  O		2
	//C	 1  2

	// Vecteur BC
	Vec3d BC = AC - AB;

	// Aire du triangle
	double area2 = 0.5 * norm(cross(AB, AC)); // = 2 * Area
	assert(area2 > 0.0 && "Triangle degeneré");

	// Facteur commun : 1 / (4 * Area)
	double coef = 1.0 / (4.0 * area2);

// Attention au signe !!! il y a un moins car dans toutes nos formules du cours on a par exemples (A-C).(B-A) qui correspondent à vect CA cross vect AB donc - vect AC cross vect AB
	S[0] = - dot(AC, BC) * coef; // AB
	S[1] = - dot(AB, BC) * coef; // AC
	S[2] = - dot(AB, AC) * coef; // BC

}
