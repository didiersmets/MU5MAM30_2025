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
	/*
    Given triangle ABC, the local mass matrix is defined as:

    integral over ABC of phi_i * phi_j dA

    where phi_i and phi_j are the shape functions associated with vertices A, B, and C.

    The triangle is represented by two edge vectors:
    AB = B - A
    AC = C - A

    Our basis functions are not defined on our triangle but on a reference triangle 
    with vertices at (0,0), (1,0), and (0,1).
    So when we compute the integral over the triangle we need to take into account
    the change of variables which introduces a Jacobian factor equal to
    the the ratio between the area of the physical triangle reference triangle ( = 1/2 ).

    1. Find the basis functions:
    phi(x,y) = a + b*x + c*y
    phi(A) for example:
    phi(0,0) = 1  => a = 1
    phi(1,0) = 0  => b = -1
    phi(0,1) = 0  => c = -1 
    phi(A) = 1 - x - y
    similarly with B and C

    2. we integrate phi(A)*phi(B) over the triangle 
      1/12 if i = j
      1/24 if i != j

    3. scale factor = Area(ABC) / Area(ref triangle) = 2*Area(ABC)
        Area(ABC) = 0.5 * || AB x AC ||

    4. final result:
        Mii = Area(ABC) / 6
        Mij = Area(ABC) / 12  if i != j
  */
  Vec3d  cross_product = cross(AB, AC);
  double area          = 0.5 * norm(cross_product);
  double diagonal_term = area / 6.0;
  double off_diagonal_term = area / 12.0;

  /*
   * Construct the mass matrix M in row-major order:
      Recall that the mass matrix is symmetric, so M[i][j] = M[j][i]
  */

  M[0] = diagonal_term;

  M[1] = off_diagonal_term;

}
