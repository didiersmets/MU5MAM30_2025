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
	Area=0.5*norm(cross(AB,AC));
	Jacobian=2*Area /* lo Jacobiano è (e i triangoli sono piccoli Area triangolo/area tr ref (0.5))*/
	/* Your implementation goes here */
	M[0]=Jacobian/12; /*(phi_i)^2 */
	M[1]=Jacobian/24; /*phi_i*phi_j  */
	
}
