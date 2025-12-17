#include <cmath>

#include "mesh/sphere.h"

int load_sphere(Mesh &m, size_t subdiv){
	load_icosahedron(m);
	return 0;
}

int load_icosahedron(Mesh &m){
	// https://en.wikipedia.org/wiki/Regular_icosahedron#Construction
	float phi = (1.0+sqrt(5.0))/2.0;
	for (int s1=-1; s1<=1; s1 +=2){
		for (int s2=-1; s2<=1; s2 +=2){
			m.positions.push_back(Vec3(0,s1 * 1, s2 * phi));
		}
	}
	for (int s1=-1; s1<=1; s1 +=2){
		for (int s2=-1; s2<=1; s2 +=2){
			m.positions.push_back(Vec3(s1 * 1,s2 * phi,0));
		}
	}
	for (int s1=-1; s1<=1; s1 +=2){
		for (int s2=-1; s2<=1; s2 +=2){
			m.positions.push_back(Vec3(s1 * phi,0,s2 * 1));
		}
	}
	m.indices.push_back(0);
	m.indices.push_back(1);
	m.indices.push_back(4);

	return 0;
}
