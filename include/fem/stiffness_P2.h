//
// Created by aurel on 20/02/2026.
//

#ifndef MU5MAM30_STIFFNESS_P2_H
#define MU5MAM30_STIFFNESS_P2_H

#endif //MU5MAM30_STIFFNESS_P2_H

#pragma once

#include <cmath>
#include "vec3.h"

/* Given a triangle ABC, computes the (symmetric) 6x6 stiffness matrix S for P2 */
inline void stiffness_P2(const Vec3d &AB, const Vec3d &AC, double S[6][6])
{
    // Produit vectoriel AB ^ AC
    double nx = AB.y * AC.z - AB.z * AC.y;
    double ny = AB.z * AC.x - AB.x * AC.z;
    double nz = AB.x * AC.y - AB.y * AC.x;

    double detJ = std::sqrt(nx*nx + ny*ny + nz*nz);

    // Calcul des coefficients c_ij = det(J) * g_ij
    double c11 =  dot(AC, AC) / detJ;
    double c22 =  dot(AB, AB) / detJ;
    double c12 = -dot(AB, AC) / detJ;

    const double S_xx[6][6] = {
        {  1./2.,  1./6.,   0., -2./3.,   0.,     0.   },
        {  1./6.,  1./2.,   0., -2./3.,   0.,     0.   },
        {   0.,     0.,     0.,   0.,     0.,     0.   },
        { -2./3., -2./3.,   0.,  4./3.,   0.,     0.   },
        {   0.,     0.,     0.,   0.,    4./3., -4./3. },
        {   0.,     0.,     0.,   0.,   -4./3.,  4./3. }
    };

    const double S_yy[6][6] = {
        {  1./2.,   0.,   1./6.,   0.,     0.,  -2./3. },
        {   0.,     0.,     0.,    0.,     0.,     0.  },
        {  1./6.,   0.,   1./2.,   0.,     0.,  -2./3. },
        {   0.,     0.,     0.,   4./3., -4./3.,   0.  },
        {   0.,     0.,     0.,  -4./3.,  4./3.,   0.  },
        { -2./3.,   0.,  -2./3.,   0.,     0.,   4./3. }
    };

    const double S_xy[6][6] = {
        {  1./2.,   0.,   1./6.,   0.,     0.,  -2./3. },
        {  1./6.,   0.,  -1./6., -2./3.,  2./3.,   0.  },
        {   0.,     0.,     0.,    0.,     0.,     0.  },
        { -2./3.,   0.,     0.,   2./3., -2./3.,  2./3.},
        {   0.,     0.,   2./3., -2./3.,  2./3., -2./3.},
        {   0.,     0.,  -2./3.,  2./3., -2./3.,  2./3.}
    };

    // Remplissage de la matrice locale S
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double s_xy_sym = S_xy[i][j] + S_xy[j][i];
            S[i][j] = c11 * S_xx[i][j] + c12 * s_xy_sym + c22 * S_yy[i][j];
        }
    }
}