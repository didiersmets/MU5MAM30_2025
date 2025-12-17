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
 * Fonctions de base P1 :
 * ---------------------
 * Les fonctions de base φ_i sont affines sur le triangle ABC :
 *     Phi_i(x,y) = a_i x + b_i y + c_i.
 * Elles sont définies par :
 *     Phi_i(S_j) = δ_ij -> φ_A(A)=1 et φ_A(B)=φ_A(C)=0 (idem pour B et C).
 *
 * Matrice de masse locale :
 * ------------------------
 * La matrice de masse M est définie par :
 *     M_ij = ∫_ABC φ_i φ_j dABC
 * Explicitement :
 *     M = [ ∫ Phi_A^2   ∫ Phi_A Phi_B   ∫ Phi_A Phi_C
 *           ∫ Phi_B Phi_A ∫ Phi_B^2     ∫ Phi_B Phi_C
 *           ∫ Phi_C Phi_A ∫ Phi_C Phi_B   ∫ Phi_C^2 ]
 *
 * Calcul pratique :
 * -----------------
 * Le triangle ABC (éventuellement plongé dans R^3) est obtenu par
 * une application affine F du triangle de référence T̂ ⊂ R^2.
 * Les fonctions de base physiques Phi_i sont les images des fonctions
 * de base de référence Phi^_i par F
 *
 * Par changement de variables :
 *     M_ij = ∫_T̂ φ̂_i φ̂_j |det(J_F)| dT̂,
 * où J_F est la matrice jacobienne de F.
 *
 * Pour un élément P1, |det(J_F)| est constant et vaut 2|ABC|,
 * ce qui conduit à la formule :
 *
 *     M = (|ABC| / 24) * [ 2  1  1
 *                         1  2  1
 *                         1  1  2 ]
 */

void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M) {
	M[0] = norm(cross(AB, AC)) / 12;
	M[1] = M[0] / 2;
}