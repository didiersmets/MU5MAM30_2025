#pragma once

#include "vec3.h"
#include "sys_utils.h"

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


 // my idea:
 // The cross product AB x AC outputs a vector that is:
 // orthogonal to both AB and AC and has magnitude equal to the area of the parallelogram
 // formed by AB and AC.
 // the magnitude is given by |AB x AC| = |AB| * |AC| * sin(theta)
 // where theta is the angle between AB and AC.
 // and this is exactly the area of the parallelogram formed by AB and AC.
 // Since the area of triangle ABC is half the area of the parallelogram 
 // we have Area_ABC = 0.5 * |AB x AC|
void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    
    // safety checks for best practices 
	// asserts are compiled out (no overhead) in release mode
	ASSERT(M != nullptr);

    ASSERT(std::isfinite(AB.x) && std::isfinite(AB.y) && std::isfinite(AB.z));
	ASSERT(std::isfinite(AC.x) && std::isfinite(AC.y) && std::isfinite(AC.z));

	// Compute area of triangle ABC
	// Area = 0.5 * |AB x AC|
	Vec3d crossProduct_AB_AC = cross(AB, AC);
	double area = 0.5 * norm(crossProduct_AB_AC);
	
	constexpr double epsilon = 1e-14;
	ASSERT(std::isfinite(area) && area > epsilon);

	// Mass matrix coefficients
	double diag_entries = area / 6.0;
	double offdiag_entries = area / 12.0;

	// Fill mass matrix M
	// the mass matrix is symmetric so Mij = Mji
	// given the dimension of M direct assignment of the values
	// is more efficient than using loops in this case
	M[0]= M[4]= M[8]= diag_entries;
	M[1] = M[2] = M[3] = M[5] = M[6] = M[7] = offdiag_entries;	

}
