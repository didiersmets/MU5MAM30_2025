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

/*

  BC = BA+AC = -AB + AC
  Area4 = Area*4 = 2 * |AB x AC|
  
  S[0] = S_00  A<->A : |BC|^2 / Area4
  S[4] = S_11  B<->B : |AC|^2 / Area4
  S[8] = S_22  C<->C : |AB|^2 / Area4
  
  S[1]=S[3] = S_01  A<->B : - (AC . BC) / Area4
  S[5]=S[4] = S_12  B<->C : - (AC . AB) / Area4
  S[2]=S[6] = S_20  C<->A :   (AB . BC) / Area4

*/



void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
  double Area4 = 2 * norm(cross(AB, AC));
  Vec3d BC = AC - AB;

  S[0] = norm2(BC) / Area4;
  S[1] = (-1) * dot(AC,BC) / Area4;
  S[2] =        dot(AB,BC) / Area4;
  
  S[4] = norm2(AC) / Area4;
  S[5] = (-1) * dot(AC,AB) / Area4;
  S[8] = norm2(AB) / Area4;
}
