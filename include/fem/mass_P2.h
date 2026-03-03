#pragma once

#include "vec3.h"

void inline mass_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
  const double area = 0.5 * norm(cross(AB, AC));
  const double c = area / 180;

  /*
    Symmetric storage: i>=j
    Index numbering:
    0: A
    1: B
    2: C
    3: AB
    4: AC
    5: BC
  */

  /* Row i=0 */
  M[0]  = 6;

  /* Row i=1 */
  M[1]  = -1;
  M[2]  = 6;

  /* Row i=2 */
  M[3]  = -1;
  M[4]  = -1;
  M[5]  = 6;

  /* Row i=3 */
  M[6]  = 0;
  M[7]  = 0;
  M[8]  = -4;
  M[9]  = 32;

  /* Row i=4 */
  M[10] = 0;
  M[11] = -4;
  M[12] = 0;
  M[13] = 16;
  M[14] = 32;

  /* Row i=5 */
  M[15] = -4;
  M[16] = 0;
  M[17] = 0;
  M[18] = 16;
  M[19] = 16;
  M[20] = 32;

  for (int i=0; i<21; ++i)
    M[i] *= c;
}
