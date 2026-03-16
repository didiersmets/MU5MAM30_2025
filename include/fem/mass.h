#include "vec3.h"
#include <cmath>

/* Given a triangle ABC, computes the (symmetric) 3x3 mass M s.t.
 *
 * M_{ij} := \int_{ABC} \phi_i \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C
 * are the shape functions of the P1 Lagrange element associated
 * to ABC.
 *
 * Idea behind computation :
 * -------------------------
 *
 * Hidden for now.
 */
void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    /* Your implementation goes here ! */

    // 1. Calculate the cross product to get the area of a triangle (2 * Area = |AB x AC|)
    double cp_x = AB[1] * AC[2] - AB[2] * AC[1];
    double cp_y = AB[2] * AC[0] - AB[0] * AC[2];
    double cp_z = AB[0] * AC[1] - AB[1] * AC[0];
    
    // 2. Calculate the area of triangle
    double area = 0.5 * std::sqrt(cp_x * cp_x + cp_y * cp_y + cp_z * cp_z);

    // 3. Fill in the mass matrix coefficients of the P1 element
    // M_{ii} = Area / 6
    // M_{ij} = Area / 12
    
    // Let M[0] store the diagonal entries, and M[1] store the off-diagonal entries
    // This needs to be consistent with the logic when assembling P1.cpp
    M[0] = area / 6.0;
    M[1] = area / 12.0;
}
