#pragma once

#include "vec3.h"
#include "sys_utils.h"

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


// For P1 elements, gradients of shape functions are constant over the triangle.
// The gradient of each shape function phi_i is computed as:
// nabla_phi_i = (n* edge_opposite_to_i) / ||n||^2
// where n = AB x AC is the normal vector to the triangle plane
// S[i][j]= Area x (nabla_phi_i · nabla_phi_j)
void inline stiffness(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{
	ASSERT(S!= nullptr);

	const Vec3d normal_vector = cross(AB, AC);
	const double norm_two_normal_vector = norm2(normal_vector);

	constexpr double epsilon = 1e-14;
	ASSERT_ALWAYS(norm_two_normal_vector > epsilon);

	const double area = 0.5 * std::sqrt(norm_two_normal_vector);
	
	// avoid recomputing division given that is slower than multiplication
	const double reciprocal_of_norm_two = 1.0 / norm_two_normal_vector;

    const Vec3d gradient_A = cross(normal_vector, AB - AC) * reciprocal_of_norm_two;
	const Vec3d gradient_B = cross(normal_vector, AC) * reciprocal_of_norm_two;
	const Vec3d gradient_C = cross(normal_vector, -AB) * reciprocal_of_norm_two;

	// S is stored as a flat 1D array
	S[0] = area * dot(gradient_A, gradient_A);  // Diagonal: A-A interaction
    S[1] = area * dot(gradient_B, gradient_B);  // Diagonal: B-B interaction
    S[2] = area * dot(gradient_C, gradient_C); // Diagonal: C-C interaction
    S[3] = area * dot(gradient_A, gradient_B); // Off-diagonal: A-B interaction
    S[4] = area * dot(gradient_B, gradient_C);  // Off-diagonal: B-C interaction
    S[5] = area * dot(gradient_C, gradient_A);  // Off-diagonal: C-A interaction
}
