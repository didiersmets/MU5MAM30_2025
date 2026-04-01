//
// Created by aurel on 20/02/2026.
//

#ifndef MU5MAM30_MASS_P2_H
#define MU5MAM30_MASS_P2_H

#endif //MU5MAM30_MASS_P2_H
#include <cmath>
#include "vec3.h"

/* Given a triangle ABC, computes the (symmetric) 6x6 mass matrix M for P2 */
inline void mass_P2(const Vec3d &AB, const Vec3d &AC, double M[6][6])
{
    // Produit vectoriel AB ^ AC pour obtenir l'aire
    double nx = AB.y * AC.z - AB.z * AC.y;
    double ny = AB.z * AC.x - AB.x * AC.z;
    double nz = AB.x * AC.y - AB.y * AC.x;

    double detJ = std::sqrt(nx*nx + ny*ny + nz*nz); // detJ = 2 * Aire
    double Aire = detJ / 2.0;

    double C = Aire / 90.0;

    // Constantes de la matrice de masse de référence
    const double M_ref[6][6] = {
        {  6.0, -1.0, -1.0,  0.0, -4.0,  0.0 },
        { -1.0,  6.0, -1.0,  0.0,  0.0, -4.0 },
        { -1.0, -1.0,  6.0, -4.0,  0.0,  0.0 },
        {  0.0,  0.0, -4.0, 32.0, 16.0, 16.0 },
        { -4.0,  0.0,  0.0, 16.0, 32.0, 16.0 },
        {  0.0, -4.0,  0.0, 16.0, 16.0, 32.0 }
    };

    // Remplissage de la matrice locale M
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            M[i][j] = M_ref[i][j] * C;
        }
    }
}