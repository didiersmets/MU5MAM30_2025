#include "mass_P2.h"
#include "stiffness_P2.h"
#include "vec3.h"
#include <cmath>

// -----------------------------------------------------------------------------
// Funciones de forma P2 y gradientes en baricéntricas
// -----------------------------------------------------------------------------

static void P2_shape_functions(double l1, double l2, double l3, double phi[6])
{
    phi[0] = l1*(2*l1 - 1);
    phi[1] = l2*(2*l2 - 1);
    phi[2] = l3*(2*l3 - 1);
    phi[3] = 4*l1*l2;
    phi[4] = 4*l2*l3;
    phi[5] = 4*l3*l1;
}

static void P2_shape_gradients(
    double l1, double l2, double l3,
    const double g1[2], const double g2[2], const double g3[2],
    double dphi[6][2])
{
    dphi[0][0] = (4*l1 - 1)*g1[0];
    dphi[0][1] = (4*l1 - 1)*g1[1];

    dphi[1][0] = (4*l2 - 1)*g2[0];
    dphi[1][1] = (4*l2 - 1)*g2[1];

    dphi[2][0] = (4*l3 - 1)*g3[0];
    dphi[2][1] = (4*l3 - 1)*g3[1];

    dphi[3][0] = 4*(l1*g2[0] + l2*g1[0]);
    dphi[3][1] = 4*(l1*g2[1] + l2*g1[1]);

    dphi[4][0] = 4*(l2*g3[0] + l3*g2[0]);
    dphi[4][1] = 4*(l2*g3[1] + l3*g2[1]);

    dphi[5][0] = 4*(l3*g1[0] + l1*g3[0]);
    dphi[5][1] = 4*(l3*g1[1] + l1*g3[1]);
}

// -----------------------------------------------------------------------------
// Cuadratura Dunavant orden 4 (6 puntos)
// -----------------------------------------------------------------------------

static const int Q = 6;
static const double lambda_Q[6][3] = {
    {0.816847572980459, 0.091576213509771, 0.091576213509771},
    {0.091576213509771, 0.816847572980459, 0.091576213509771},
    {0.091576213509771, 0.091576213509771, 0.816847572980459},
    {0.108103018168070, 0.445948490915965, 0.445948490915965},
    {0.445948490915965, 0.108103018168070, 0.445948490915965},
    {0.445948490915965, 0.445948490915965, 0.108103018168070}
};

static const double w_Q[6] = {
    0.109951743655322,
    0.109951743655322,
    0.109951743655322,
    0.223381589678011,
    0.223381589678011,
    0.223381589678011
};

// -----------------------------------------------------------------------------
// MASS MATRIX P2
// -----------------------------------------------------------------------------

void mass_P2_custom(const Vec3d &AB, const Vec3d &AC, double Mloc[36])
{
    // Coordenadas del triángulo
    double x1 = 0, y1 = 0;
    double x2 = AB.x, y2 = AB.y;
    double x3 = AC.x, y3 = AC.y;

    double detT = (x2*y3 - x3*y2);
    double area = 0.5 * fabs(detT);

    // Gradientes de baricéntricas
    double g1[2] = { (y2 - y3)/detT, (x3 - x2)/detT };
    double g2[2] = { (y3 - y1)/detT, (x1 - x3)/detT };
    double g3[2] = { (y1 - y2)/detT, (x2 - x1)/detT };

    for (int i = 0; i < 36; ++i)
        Mloc[i] = 0.0;

    for (int q = 0; q < Q; ++q) {
        double l1 = lambda_Q[q][0];
        double l2 = lambda_Q[q][1];
        double l3 = lambda_Q[q][2];

        double phi[6];
        P2_shape_functions(l1, l2, l3, phi);

        double w = w_Q[q] * area;

        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                Mloc[i*6 + j] += phi[i] * phi[j] * w;
    }
}

// -----------------------------------------------------------------------------
// STIFFNESS MATRIX P2
// -----------------------------------------------------------------------------

void stiffness_P2_custom(const Vec3d &AB, const Vec3d &AC, double Sloc[36])
{
    double x1 = 0, y1 = 0;
    double x2 = AB.x, y2 = AB.y;
    double x3 = AC.x, y3 = AC.y;

    double detT = (x2*y3 - x3*y2);
    double area = 0.5 * fabs(detT);

    double g1[2] = { (y2 - y3)/detT, (x3 - x2)/detT };
    double g2[2] = { (y3 - y1)/detT, (x1 - x3)/detT };
    double g3[2] = { (y1 - y2)/detT, (x2 - x1)/detT };

    for (int i = 0; i < 36; ++i)
        Sloc[i] = 0.0;

    for (int q = 0; q < Q; ++q) {
        double l1 = lambda_Q[q][0];
        double l2 = lambda_Q[q][1];
        double l3 = lambda_Q[q][2];

        double dphi[6][2];
        P2_shape_gradients(l1, l2, l3, g1, g2, g3, dphi);

        double w = w_Q[q] * area;

        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                Sloc[i*6 + j] += (dphi[i][0]*dphi[j][0] +
                                  dphi[i][1]*dphi[j][1]) * w;
    }
}
