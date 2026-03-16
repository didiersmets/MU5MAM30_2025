#pragma once

#include "vec3.h"
#include <cmath>

/* Given a triangle ABC, computes the (symmetric) 3x3 stiffness matrix S s.t.
 *
 * S_{ij} := \int_{ABC} \nabla \phi_i \cdot \nabla \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C
 * are the shape functions of the P1 Lagrange element associated to ABC.
 *
 * Output: the six coefficients S_{00} S_{11} S_{22} S_{01} S_{12} S_{20},
 * corresponding to the interactions A<->A, B<->B, C<->C, A<->B, B<->C, C<->A
 */
// void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
// {
// 	/* Your implementation goes here */
// }


void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
    // 1. Define the three edge vectors (for geometric computation of gradient inner products)
    // the 3 numbers in {} are the 3 dimentions
    // Edge vector corresponding to vertex A is BC = AC - AB
    Vec3d v0 = {AC[0] - AB[0], AC[1] - AB[1], AC[2] - AB[2]}; 
    // Edge vector corresponding to vertex B is CA = -AC
    Vec3d v1 = {-AC[0], -AC[1], -AC[2]}; 
    // Edge vector corresponding to vertex C is AB = AB
    Vec3d v2 = {AB[0], AB[1], AB[2]};

    // 2. Compute cross product to obtain the triangle area
    double cp_x = AB[1] * AC[2] - AB[2] * AC[1];
    double cp_y = AB[2] * AC[0] - AB[0] * AC[2];
    double cp_z = AB[0] * AC[1] - AB[1] * AC[0];
    
    // area2 = 2 * Area (area of the parallelogram)
    double area2 = std::sqrt(cp_x * cp_x + cp_y * cp_y + cp_z * cp_z);

    // Precompute coefficient 1 / (2 * area2), equivalent to 1 / (4 * Area)
    double inv_area4 = 1.0 / (2.0 * area2);

    // 3. Fill diagonal terms S_ii = (L_i · L_i) / (4 * Area)
    // S00: vertex A self-interaction
    S[0] = (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]) * inv_area4;
    // S11: vertex B self-interaction
    S[1] = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) * inv_area4;
    // S22: vertex C self-interaction
    S[2] = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]) * inv_area4;

    // 4. Fill off-diagonal terms S_ij = (L_i · L_j) / (4 * Area)
    // S01: interaction between A and B
    S[3] = (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]) * inv_area4;
    // S12: interaction between B and C
    S[4] = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) * inv_area4;
    // S20: interaction between C and A
    S[5] = (v2[0] * v0[0] + v2[1] * v0[1] + v2[2] * v0[2]) * inv_area4;
}