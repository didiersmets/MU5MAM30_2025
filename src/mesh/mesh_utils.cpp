#include "mesh_utils.h"

double area_triangle(const Vec3 &AB,const Vec3 &AC){
	double area = 0.5*(AB.x * AC.y - AB.y * AC.x);
	return area;
}
