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

void inline stiffness_v2(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
    // Étape 1 : Calcul des produits scalaires et normes carrées
    const double AB_sq = norm2(AB);      // |AB|²
    const double AC_sq = norm2(AC);      // |AC|²
    const double AB_dot_AC = dot(AB, AC); // AB·AC
    
    // Étape 2 : Calcul de l'aire du triangle via le produit vectoriel
    // Aire² = (|AB|² × |AC|² - (AB·AC)²) / 4
    const double four_area_sq = AB_sq * AC_sq - AB_dot_AC * AB_dot_AC;
    const double inv_four_area = 1.0 / sqrt(four_area_sq);
    
    // Étape 3 : Calcul du vecteur BC = AC - AB
    const double BC_sq = AB_sq + AC_sq - 2.0 * AB_dot_AC; // |BC|² = |AC-AB|²
    
    // Étape 4 : Coefficients diagonaux S_ii = |côté_opposé|² / (4×Aire)
    S[0] = BC_sq * inv_four_area;  // S_AA : côté opposé = BC
    S[1] = AC_sq * inv_four_area;  // S_BB : côté opposé = AC
    S[2] = AB_sq * inv_four_area;  // S_CC : côté opposé = AB
    
    // Étape 5 : Coefficients hors-diagonal S_ij = -cot(angle) / 2
    // où angle est l'angle au sommet opposé à l'arête (i,j)
    S[3] = -AB_dot_AC * inv_four_area;              // S_AB : angle en C
    S[4] = -(AC_sq - AB_dot_AC) * inv_four_area;    // S_BC : angle en A
    S[5] = -(AB_sq - AB_dot_AC) * inv_four_area;    // S_CA : angle en B
}

/* Alternative avec structure explicite pour plus de clarté */
struct StiffnessMatrix {
    double AA, BB, CC;  // Termes diagonaux
    double AB, BC, CA;  // Termes hors-diagonal
};

inline StiffnessMatrix stiffness_structured(const Vec3d &AB, const Vec3d &AC)
{
    const double AB_sq = norm2(AB);
    const double AC_sq = norm2(AC);
    const double AB_dot_AC = dot(AB, AC);
    
    // Aire via formule de Héron vectorielle
    const double inv_four_area = 1.0 / sqrt(AB_sq * AC_sq - AB_dot_AC * AB_dot_AC);
    
    // BC² = |AC - AB|²
    const double BC_sq = AB_sq + AC_sq - 2.0 * AB_dot_AC;
    
    StiffnessMatrix S;
    
    // Diagonale : formule du côté opposé
    S.AA = BC_sq * inv_four_area;
    S.BB = AC_sq * inv_four_area;
    S.CC = AB_sq * inv_four_area;
    
    // Hors-diagonal : formule de la cotangente
    S.AB = -AB_dot_AC * inv_four_area;
    S.BC = -(AC_sq - AB_dot_AC) * inv_four_area;
    S.CA = -(AB_sq - AB_dot_AC) * inv_four_area;
    
    return S;
}
	
