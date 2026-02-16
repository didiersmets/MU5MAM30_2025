#pragma once

#include "vec3.h"

// Matriz de rigidez P2 (6x6) en un triángulo ABC
// Entrada:
//   AB = B - A
//   AC = C - A
// Salida:
//   S[36] = matriz 6x6 en orden fila-major: S[i*6 + j]
inline void stiffness_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
    // 1. Área del triángulo
    double area2 = norm(cross(AB, AC));   // = 2 * area
    double area  = 0.5 * area2;

    // 2. Matriz de rigidez P2 en el triángulo de referencia
    static const double Kref[36] = {
         3, -1, -1, -4,  2,  2,
        -1,  3, -1,  2, -4,  2,
        -1, -1,  3,  2,  2, -4,
        -4,  2,  2, 16, -8, -8,
         2, -4,  2, -8, 16, -8,
         2,  2, -4, -8, -8, 16
    };

    // 3. Escalado correcto: 1 / (24 * area)
    double scale = 1.0 / (24.0 * area);

    // 4. Construir la matriz final
    for (int i = 0; i < 36; ++i)
        S[i] = scale * Kref[i];
}
