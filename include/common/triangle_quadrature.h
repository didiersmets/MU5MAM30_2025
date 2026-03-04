#pragma once

#include "vec3.h"

template<typename F>
double integrate_triangle_deg5(const Vec3 &A,
                               const Vec3 &B,
                               const Vec3 &C,
                               F integrand)
{
  /*
    Integrate over the triangle ABC using the degree-5 exact 7-point Dunavant
    quadrature formula
  */

  Vec3 AB = B - A;
  Vec3 AC = C - A;
  Vec3 N  = cross(AB,AC);
  double area = 0.5 * norm(N);

  const double w0 = 0.225;
  const double w1 = 0.1323941527;
  const double w2 = 0.1259391805;

  const double a1 = 0.0597158717;
  const double b1 = 0.4701420641;

  const double a2 = 0.7974269853;
  const double b2 = 0.1012865073;

  double result = 0.0;

  result += w0 * integrand(1.0/3, 1.0/3, 1.0/3);

  result += w1 * integrand(a1,b1,b1);
  result += w1 * integrand(b1,a1,b1);
  result += w1 * integrand(b1,b1,a1);

  result += w2 * integrand(a2,b2,b2);
  result += w2 * integrand(b2,a2,b2);
  result += w2 * integrand(b2,b2,a2);

  return area * result;
}
