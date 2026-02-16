#pragma once

#include "vec3.h"

// Matriz de masa P2 (6x6) en un triángulo ABC
// Entrada:
//   AB = B - A
//   AC = C - A
// Salida:
//   M[36] = matriz 6x6 en orden fila-major: M[i*6 + j]
inline void mass_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict M)
{
    // 1. Área del triángulo
    double area2 = norm(cross(AB, AC));   // = 2 * area
    double area  = 0.5 * area2;

    // 2. Matriz de masa P2 en el triángulo de referencia (sin área)
    //    M_ref = (1/180) * [...]
    static const double Mref[36] = {
        // fila 0
         6, -1, -1,  0,  0,  0,
        // fila 1
        -1,  6, -1,  0,  0,  0,
        // fila 2
        -1, -1,  6,  0,  0,  0,
        // fila 3
         0,  0,  0, 32, 16, 16,
        // fila 4
         0,  0,  0, 16, 32, 16,
        // fila 5
         0,  0,  0, 16, 16, 32
    };

    // 3. Escalar por el factor correcto: area / 180
    double scale = area / 180.0;

    // 4. Construir la matriz final
    for (int i = 0; i < 36; ++i)
        M[i] = scale * Mref[i];
}
