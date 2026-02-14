#pragma once

#include "vec3.h"

static void build_L(const Vec3d &AB, const Vec3d &AC, double *__restrict L)
{
	Vec3d CB = AB - AC;

	/* Computation of AB x AC */
	Vec3d n = cross(AB, AC);
	double K2 = norm2(n);

	L[0] = norm2(CB) / K2; // A <-> A
	L[1] = norm2(AC) / K2; // A <-> A
	L[2] = norm2(AB) / K2; // A <-> A

	L[3] = dot(CB, AC) / K2;  // A <-> B
	L[4] = -dot(CB, AB) / K2; // A <-> C
	L[5] = -dot(AC, AB) / K2; // B <-> C
}

void inline stiffness_P2(const Vec3d &AB, const Vec3d &AC, double *__restrict S)
{

	double L[6];
	build_L(AB, AC, L);

	/* Computation of ||AB x AC|| */
	double K_6 = norm(cross(AB, AC)) / 6;

	if (K_6 == 0)
		throw std::runtime_error("The two vectors must be linearly independent.");
	else
	{
		/* Computation of the coefficients in the bloc vertex-vertex */
		S[0] = K_6 * 3 * L[0];
		S[1] = -K_6 * L[3];
		S[2] = -K_6 * L[4];

		S[7] = K_6 * 3 * L[1];
		S[8] = -K_6 * L[5];

		S[14] = K_6 * 3 * L[2];

		/* Computation of the coefficients in the bloc vertex-edge */
		S[3] = K_6 * 4 * L[3];
		S[4] = 0;
		S[5] = K_6 * 4 * L[4];

		S[9] = K_6 * 4 * L[3];
		S[10] = K_6 * 4 * L[5];
		S[11] = 0;

		S[15] = 0;
		S[16] = K_6 * 4 * L[5];
		S[17] = K_6 * 4 * L[4];

		/* Computation of the coefficients in the bloc edge-edge */
		S[21] = K_6 * 8 * (L[0] - L[5]);
		S[22] = K_6 * 8 * L[4];
		S[23] = K_6 * 8 * L[5];

		S[28] = K_6 * 8 * (L[1] - L[4]);
		S[29] = K_6 * 8 * L[3];

		S[35] = K_6 * 8 * (L[2] - L[3]);

		/* Computation of the symmetric part */
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < i; j++)
				S[6 * i + j] = S[6 * j + i];
	}
}
