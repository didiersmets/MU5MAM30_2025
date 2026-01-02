
#pragma once

#include "common/vec3.h"

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
 * Idea behind computation : see class notes
 *
 */
void inline stiffness(const Vec3d& AB, const Vec3d& AC, double* __restrict S)
{
  /* Computation of ||AB x AC|| */
  double det = norm(cross(AB, AC));

  Vec3d BC = AC - AB;

  // the diagonal part
  S[0] = norm2(BC) / (2 * det);  // A <-> A
  S[1] = norm2(AC) / (2 * det);  // B <-> B
  S[2] = norm2(AB) / (2 * det);  // C <-> C

  // the off-diagonal part
  S[3] = -dot(AC, BC) / (2 * det);  // A <-> B
  S[4] = -dot(AB, AC) / (2 * det);  // B <-> C
  S[5] = dot(AB, BC) / (2 * det);   // C <-> A
}