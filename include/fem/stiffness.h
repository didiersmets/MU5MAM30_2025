#pragma once

#include "../common/vec3.h"

#include <iostream>
#include "mass.h"      

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


void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
    // 1) Área del triángulo
    const double area = 0.5 * norm(cross(AB, AC));

    // 2) Construir Jacobiano J = [AB AC]
    // En 2D: J = |AB_x AC_x|
    //           |AB_y AC_y|
    double J[2][2] = {
        {AB[0], AC[0]},
        {AB[1], AC[1]}
    };

    // 3) Calcular inversa transpuesta de J (J^-T)
    const double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    double JTinv[2][2] = {
        {  J[1][1]/detJ, -J[1][0]/detJ},
        {-J[0][1]/detJ,  J[0][0]/detJ}
    };

    // 4) Gradientes en coordenadas baricéntricas
    double gradAlphaBeta[3][2] = {
        { 1.0,  0.0},  // grad φ_A
        { 0.0,  1.0},  // grad φ_B
        {-1.0, -1.0}   // grad φ_C
    };

    // 5) Transformar a gradientes en coordenadas físicas
    Vec3d grad[3]; // guardamos como 3D pero solo x,y se usan
    for(int i=0; i<3; ++i)
    {
        grad[i][0] = JTinv[0][0]*gradAlphaBeta[i][0] + JTinv[0][1]*gradAlphaBeta[i][1]; // x
        grad[i][1] = JTinv[1][0]*gradAlphaBeta[i][0] + JTinv[1][1]*gradAlphaBeta[i][1]; // y
        grad[i][2] = 0.0; // z = 0
    }

    // 6) Calcular los productos punto y multiplicar por el área
    const double S_AA = (grad[0][0]*grad[0][0] + grad[0][1]*grad[0][1]) * area;
    const double S_BB = (grad[1][0]*grad[1][0] + grad[1][1]*grad[1][1]) * area;
    const double S_CC = (grad[2][0]*grad[2][0] + grad[2][1]*grad[2][1]) * area;

    const double S_AB = (grad[0][0]*grad[1][0] + grad[0][1]*grad[1][1]) * area;
    const double S_BC = (grad[1][0]*grad[2][0] + grad[1][1]*grad[2][1]) * area;
    const double S_CA = (grad[2][0]*grad[0][0] + grad[2][1]*grad[0][1]) * area;

    // 7) Guardar en arreglo según convención
    S[0] = S_AA;
    S[1] = S_BB;
    S[2] = S_CC;
    S[3] = S_AB;
    S[4] = S_BC;
    S[5] = S_CA;
}
