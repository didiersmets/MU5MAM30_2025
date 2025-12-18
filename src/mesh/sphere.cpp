#include <cmath>

#include "sphere.h"

int load_sphere(Mesh &m, size_t subdiv){
	load_icosahedron(m);
	return 0;
}

void add_face(Mesh &m, uint32_t f1, uint32_t f2, uint32_t f3){
	m.indices.push_back(f1);
	m.indices.push_back(f2);
	m.indices.push_back(f3);
}

int load_icosahedron(Mesh &m){
	// https://en.wikipedia.org/wiki/Regular_icosahedron#Construction
	float phi = (1.0+sqrt(5.0))/2.0;
	m.positions.push_back(Vec3(-1,-phi,0));
	m.positions.push_back(Vec3(0,-1,phi));
	m.positions.push_back(Vec3(-phi,0,1));
	m.positions.push_back(Vec3(-phi,0,-1));
	m.positions.push_back(Vec3(0,-1,-phi));
	m.positions.push_back(Vec3(1,-phi,0));
	m.positions.push_back(Vec3(0,1,phi));
	m.positions.push_back(Vec3(-1,phi,0));
	m.positions.push_back(Vec3(0,1,-phi));
	m.positions.push_back(Vec3(phi,0,-1));
	m.positions.push_back(Vec3(phi,0,1));
	m.positions.push_back(Vec3(1,phi,0));
	add_face(m,1,2,3);
	add_face(m,2,1,6);
	add_face(m,1,3,4);
	add_face(m,1,4,5);
	add_face(m,1,5,6);
	add_face(m,2,6,11);
	add_face(m,3,2,7);
	add_face(m,4,3,8);
	add_face(m,5,4,9);
	add_face(m,6,5,10);
	add_face(m,2,11,7);
	add_face(m,3,7,8);
	add_face(m,4,8,9);
	add_face(m,5,9,10);
	add_face(m,6,10,11);
	add_face(m,7,11,12);
	add_face(m,8,7,12);
	add_face(m,9,8,12);
	add_face(m,10,9,12);
	add_face(m,11,10,12);
	return 0;
}
