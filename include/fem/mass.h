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
	/* Your implementation goes here ! */
		  // Compute the area of the triangle ABC.
  Vec3d  cross_product = cross(AB, AC);
  double area          = 0.5 * norm(cross_product);

  // Mii = Area / 6.
  double diagonal_term = area / 6.0;

  // Mij = Area / 12.
  double off_diagonal_term = area / 12.0;

  /*
   * Construct the mass matrix M in row-major order:
      Recall that the mass matrix is symmetric, so M[i][j] = M[j][i]
  */

  M[0] = diagonal_term;

  M[1] = off_diagonal_term;

}
