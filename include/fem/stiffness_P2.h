#pragma once

#include "vec3.h"

void inline stiffness_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
  const double area = 0.5 * norm(cross(AB, AC));
  const double c = 1 / (12 * area);

  const Vec3d BC = AC - AB;

  const double G11 = norm2(BC);
  const double G22 = norm2(AC);
  const double G33 = norm2(AB);
  const double G12 = dot(BC, -AC);
  const double G13 = dot(BC, AB);
  const double G23 = dot(-AC, AB);

  /* Row i=0 */
  S[0]  = 3*G11;

  /* Row i=1 */
  S[1]  = -G12;
  S[2]  = 3*G22;

  /* Row i=2 */
  S[3]  = -G13;
  S[4]  = -G23;
  S[5]  = 3*G33;

  /* Row i=3 */
  S[6]  = 4*G12;
  S[7]  = 4*G12;
  S[8]  = 0;
  S[9]  = 8*(G11-G23);

  /* Row i=4 */
  S[10] = 4*G13;
  S[11] = 0;
  S[12] = 4*G13;
  S[13] = 8*G23;
  S[14] = 8*(G33-G12);

  /* Row i=5 */
  S[15] = 0;
  S[16] = 4*G23;
  S[17] = 4*G23;
  S[18] = 8*G13;
  S[19] = 8*G12;
  S[20] = 8*(G22-G13);

  for (int i=0; i<21; ++i)
    S[i] *= c;
}
