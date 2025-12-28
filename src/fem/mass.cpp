#include "mass.h"
#include "mesh_utils.h"


MassCoef mass(const Vec3 &AB, const Vec3 &AC)
{
	double area = area_triangle(AB,AC);
	return MassCoef {
		.diag = area/6,
		.offdiag = area/12
	};
}
