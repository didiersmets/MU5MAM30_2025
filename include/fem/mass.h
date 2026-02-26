#pragma once

#include "vec3.h"

/* Given a triangle ABC, computes the (symmetric) 3x3 mass M s.t.
 *
 *   M_{ij} := \int_{ABC} \phi_i \phi_j
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C
 * are the shape functions of the P1 Lagrange element associated
 * to ABC.
 *
 * Idea behind computation :
 * -------------------------
 *
 * Hiden for now.
 */
void inline mass(const Vec3d &AB, const Vec3d &AC, double *__restrict M) {
	double Area = 0.5 * std::sqrt(std::pow(norm(AB), 2)*std::pow(norm(AC), 2)-std::pow(dot(AB, AC), 2));

	for (int i= 0; i < 3; i++) {
		for (int j= 0; j < 3; j++) {
			if (i != j) {
				M[3*i + j] = Area/12;
			}
			else {
				M[3*i + j] = Area/6;
			}
		}
	}
}
