#include <stdint.h>
#include <stdio.h>
#include <cmath>
#include <utility>

#include "cube.h"
#include "mesh.h"
#include "vec3.h"

#include "sphere.h"

void push_triangle(Mesh &m, uint32_t f1, uint32_t f2, uint32_t f3){
	m.indices.push_back(f1);
	m.indices.push_back(f2);
	m.indices.push_back(f3);
}

void subdiv_step(Mesh &m){
	size_t triangle_count = m.triangle_count();
	size_t vertex_count = m.vertex_count();
	TArray<uint32_t> previous_indices = std::move(m.indices);
	TArray<int32_t> midpoint(vertex_count*vertex_count, -1);
	for (size_t t=0;t<triangle_count;t++){
		uint32_t i[3] = {previous_indices[t*3],previous_indices[(t*3)+1],previous_indices[(t*3)+2]};
		Vec3 v[3] = {m.positions[i[0]],m.positions[i[1]],m.positions[i[2]]};
		uint32_t childs[3];
		for (size_t edge=0; edge<3; edge++){
			size_t mid_i = (i[edge]<=i[(edge+1)%3]) ? vertex_count*i[edge]+i[(edge+1)%3] : vertex_count*i[(edge+1)%3]+i[edge];
			if (midpoint[mid_i] == -1){
				// construct new child
				Vec3 child = normalized(v[edge]+v[(edge+1)%3]);
				uint32_t childi = m.vertex_count();
				m.positions.push_back(child);
				childs[edge] = childi;
				midpoint[mid_i] = childi;
			}
			else {
				// child already constructed by other parent
				childs[edge] = midpoint[mid_i];
			}
		}
		push_triangle(m,i[0],childs[0],childs[2]);
		push_triangle(m,i[1],childs[1],childs[0]);
		push_triangle(m,i[2],childs[2],childs[1]);
		push_triangle(m,childs[0],childs[1],childs[2]);
	}
}

int load_sphere(Mesh &m, size_t subdiv){
	load_icosahedron(m);
	for (size_t i=0;i<subdiv;i++){
		subdiv_step(m);
	}
	return 0;
}

int load_tetrahedron(Mesh &m){
	// https://en.wikipedia.org/wiki/Regular_tetrahedron#Cartesian_coordinates
	float isqrt = 1.0/sqrt(2.0);
	m.positions.push_back(Vec3(1,0,-isqrt));
	m.positions.push_back(Vec3(-1,0,-isqrt));
	m.positions.push_back(Vec3(0,-1,isqrt));
	m.positions.push_back(Vec3(0,1,isqrt));
	push_triangle(m,0,1,2);
	push_triangle(m,0,1,3);
	push_triangle(m,2,3,0);
	push_triangle(m,2,3,1);
	// normalize everything to get a radius 1 sphere
	for (size_t i=0; i<m.vertex_count();i++){
		m.positions[i] = normalized(m.positions[i]);
	}
	return 0;
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
	push_triangle(m,1-1,2-1,3-1);
	push_triangle(m,2-1,1-1,6-1);
	push_triangle(m,1-1,3-1,4-1);
	push_triangle(m,1-1,4-1,5-1);
	push_triangle(m,1-1,5-1,6-1);
	push_triangle(m,2-1,6-1,11-1);
	push_triangle(m,3-1,2-1,7-1);
	push_triangle(m,4-1,3-1,8-1);
	push_triangle(m,5-1,4-1,9-1);
	push_triangle(m,6-1,5-1,10-1);
	push_triangle(m,2-1,11-1,7-1);
	push_triangle(m,3-1,7-1,8-1);
	push_triangle(m,4-1,8-1,9-1);
	push_triangle(m,5-1,9-1,10-1);
	push_triangle(m,6-1,10-1,11-1);
	push_triangle(m,7-1,11-1,12-1);
	push_triangle(m,8-1,7-1,12-1);
	push_triangle(m,9-1,8-1,12-1);
	push_triangle(m,10-1,9-1,12-1);
	push_triangle(m,11-1,10-1,12-1);
	// normalize everything to get a radius 1 sphere
	for (size_t i=0; i<m.vertex_count();i++){
		m.positions[i] = normalized(m.positions[i]);
	}
	return 0;
}

int load_spherical_cube(Mesh &m, size_t subdiv)
{
	if (int res = load_cube(m, subdiv))
		return (res);

	Vec3 *pos = m.positions.data;
	size_t vtx_count = m.positions.size;
	for (size_t i = 0; i < vtx_count; ++i) {
		pos[i] = normalized(pos[i]);
	}

	return (0);
}

int load_spherical_nested_cube(Mesh &m, size_t subdiv)
{
	if (int res = load_cube_nested_dissect(m, subdiv))
		return (res);

	Vec3 *pos = m.positions.data;
	size_t vtx_count = m.positions.size;
	for (size_t i = 0; i < vtx_count; ++i) {
		pos[i] = normalized(pos[i]);
	}

	return (0);
}
