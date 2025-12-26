#include "mass.h"
#include "mesh_utils.h"


MassCoef inline mass(const Vec3d &AB, const Vec3d &AC);
{
	double area = area_triangle(AB,AC);
	return MassCoef {
		.diag = area/6;
		.offdiag = area/12;
	};
}
