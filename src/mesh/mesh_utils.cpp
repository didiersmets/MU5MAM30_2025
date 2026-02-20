#include "mesh_utils.h"

double area_triangle(const Vec3d &AB,const Vec3d &AC){
	double area = 0.5*norm(cross(AB,AC));
	return area;
}
