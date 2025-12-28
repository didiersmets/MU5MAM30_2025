#include "stiffness.h"
#include "mesh_utils.h"

StiffnessCoef stiffness(const Vec3 &AB, const Vec3 &AC)
{
	double area = area_triangle(AB,AC);
	/* Your implementation goes here */
	double AB2 = norm2(AB);
	double AC2 = norm2(AC);
	double BC2 = norm2(AC-AB);
	return StiffnessCoef {
		.S00 = BC2/(4 * area),
		.S11 = AC2/(4 * area),
		.S22 = AB2/(4 * area),
		.S01 = dot(AB-AC,AC)/(4*area),
		.S12 = -dot(AC,AB)/(4*area),
		.S20 = dot(AB,AC-AB)/(4*area)
	};
}
