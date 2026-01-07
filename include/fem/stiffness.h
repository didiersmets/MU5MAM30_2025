#pragma once

#include "../common/vec3.h"

#include <iostream>
#include "mass.h" 
#pragma once
#include <cmath>     

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



// Utilidades básicas
static inline Vec3d cross(const Vec3d& a, const Vec3d& b) {
    return Vec3d{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static inline double norm(const Vec3d& v) {
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static inline Vec3d normalize(const Vec3d& v) {
    double n = norm(v);
    return (n > 0.0) ? Vec3d{v.x/n, v.y/n, v.z/n} : Vec3d{0.0,0.0,0.0};
}

/*
 * Given triangle ABC via AB=B-A and AC=C-A,
 * compute stiffness entries:
 *   S00 S11 S22 S01 S12 S20  (row-major pairs AA, BB, CC, AB, BC, CA)
 */
void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
    // Área 3D
    const Vec3d n = cross(AB, AC);
    const double area = 0.5 * norm(n);
    if (area <= 0.0) {
        for (int k = 0; k < 6; ++k) S[k] = 0.0;
        return;
    }

    // Base ortonormal del plano: t1, t2
    // t1 paralelo a AB, t2 en el plano, ortogonal a t1
    Vec3d t1 = normalize(AB);
    Vec3d t2 = cross(normalize(n), t1); // = (n̂ × t1), está en el plano y ortogonal a t1
    // Proyecciones 2D de AB y AC en esta base
    const double ABx = AB.x * t1.x + AB.y * t1.y + AB.z * t1.z;           // ~ ||AB||
    const double ABy = AB.x * t2.x + AB.y * t2.y + AB.z * t2.z;           // ~ 0
    const double ACx = AC.x * t1.x + AC.y * t1.y + AC.z * t1.z;
    const double ACy = AC.x * t2.x + AC.y * t2.y + AC.z * t2.z;

    // Jacobiano 2D y su determinante (equivale al doble del área firmado en la base local)
    const double detJ = ABx * ACy - ACx * ABy;
    if (std::abs(detJ) <= 0.0) { // degenerado en la base local
        for (int k = 0; k < 6; ++k) S[k] = 0.0;
        return;
    }

    // Gradientes en 2D (gA, gB, gC)
    const double invDet = 1.0 / detJ;

    // gA = J^{-T} * (-1,-1)
    const double gAx = ( -ACy + ABy ) * invDet;
    const double gAy = (  ACx - ABx ) * invDet;

    // gB = J^{-T} * (1,0)
    const double gBx = (  ACy ) * invDet;
    const double gBy = ( -ACx ) * invDet;

    // gC = J^{-T} * (0,1)
    const double gCx = ( -ABy ) * invDet;
    const double gCy = (  ABx ) * invDet;

    // Productos punto de gradientes
    const double gA_gA = gAx*gAx + gAy*gAy;
    const double gB_gB = gBx*gBx + gBy*gBy;
    const double gC_gC = gCx*gCx + gCy*gCy;

    const double gA_gB = gAx*gBx + gAy*gBy;
    const double gB_gC = gBx*gCx + gBy*gCy;
    const double gC_gA = gCx*gAx + gCy*gAy;

    // Entradas S_ij = (∇φi · ∇φj) * area
    S[0] = gA_gA * area; // S00 (A,A)
    S[1] = gB_gB * area; // S11 (B,B)
    S[2] = gC_gC * area; // S22 (C,C)

    S[3] = gA_gB * area; // S01 (A,B)
    S[4] = gB_gC * area; // S12 (B,C)
    S[5] = gC_gA * area; // S20 (C,A)
}

