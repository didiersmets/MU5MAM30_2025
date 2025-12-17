#pragma once

#include "vec3.h"

#include <stdexcept>

/**
 * Computes the symmetric 3x3 stiffness matrix S for a planar P1 finite element triangle in
 * R^3.
 *
 * The computation uses the analytic result derived from the change of coordinates (isoparametric
 * mapping) and the constant Metric Tensor G = J^T J. The final formula simplifies to geometric
 * products (dot, norm2) scaled by the triangle's Area(T).
 *
 * S_{ij} := \int_{T} \nabla \phi_i \cdot \nabla \phi_j \ d\Omega
 *
 * where \phi_0 := \phi_A, \phi_1 := \phi_B, \phi_2 := \phi_C are the P1 shape functions.
 *
 *  Vector from vertex A to B (B - A) in R^3.
 *  Vector from vertex A to C (C - A) in R^3.
 *  Output array (size 9) to store the stiffness matrix coefficients in row-major order.
 */
void inline stiffness(const Vec3d& AB, const Vec3d& AC, double* __restrict S)
{
	/* Computation of ||AB x AC|| */
	double det = norm(cross(AB, AC));
	assert(det != 0);
	/* Computation of the upper coefficients */
	S[0] = norm2(AC - AB) / (2 * det);
	S[1] = -dot(AC, AC - AB) / (2 * det);
	S[2] = dot(AB, AC - AB) / (2 * det);
	S[4] = norm2(AC) / (2 * det);
	S[5] = -dot(AB, AC) / (2 * det);
	S[8] = norm2(AB) / (2 * det);
	
	/* Computation of the lower coefficients */
	for (int j = 0; j < 3; j++)
	for (int i = j + 1; i < 3; i++)
		S[3 * i + j] = S[3 * j + i];
}