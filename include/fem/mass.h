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
	// Triangle area is half the magnitude of the cross product of its edges.
	const double area = 0.5 * norm(cross(AB, AC));

	// Consistent P1 mass matrix on a triangle:
	// [2 1 1] * area / 12. Store only the unique entries in the order
	// A-A, B-B, C-C, A-B, B-C, C-A.
	const double m_diag = area / 6.0;   // 2/12 * area
	const double m_off = area / 12.0;   // 1/12 * area

	M[0] = m_diag; // phi_A * phi_A
	M[1] = m_off;  // phi_C * phi_A
}
