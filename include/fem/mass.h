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
    // cross = AB x AC
    Vec3d cross;
    cross.x = AB.y * AC.z - AB.z * AC.y;
    cross.y = AB.z * AC.x - AB.x * AC.z;
    cross.z = AB.x * AC.y - AB.y * AC.x;

    double cross_norm = std::sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
    double area = 0.5 * cross_norm; //l'aire du triangle ABC vaut 1/2 |AB x AC|

    double m_diag = area / 6.0;   // int phi_i^2
    double m_off  = area / 12.0;  // int phi_i phi_j  (i != j)
    M[0] = m_diag; M[1] = m_off; 
}

