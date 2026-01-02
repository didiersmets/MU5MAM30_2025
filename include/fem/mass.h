#pragma once

#include "common/vec3.h"

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
 * see pdf FEM_mass_stiffness.pdf
 *
 * * The mass matrix M is given by:
 * M_ii = Area(ABC) / 6
 * M_ij = Area(ABC) / 12 for i != j
 *
 * Area(ABC) = 0.5 * ||AB x AC|| (assuming ABC are in 3D, and the area is calculated
 * in the plane where the triangle lies, based on 3D vectors)

 */
void inline mass(const Vec3d& AB, const Vec3d& AC, double* __restrict M)
{
  // Compute the area of the triangle ABC.
  Vec3d  cross_product = cross(AB, AC);
  double area          = 0.5 * norm(cross_product);

  // Mii = Area / 6.
  double diagonal_term = area / 6.0;

  // Mij = Area / 12.
  double off_diagonal_term = area / 12.0;

  // Diagonal Term
  M[0] = diagonal_term;

  // Non diagonal Term
  M[1] = off_diagonal_term;
}
