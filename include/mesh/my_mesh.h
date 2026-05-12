#pragma once

#include <iostream>
#include <string.h>
#include "vec2.h"
#include "vec3.h"


struct MyMesh {
	size_t vtx_count;
	size_t tri_count;
	struct TVec2<double> *vertices;
	struct TVec3<size_t> *triangles;
	size_t *boundary;
	size_t boundary_count;
};

void build_square_vertices(struct TVec2<double> *vert, size_t N, double a, double b);
void build_square_triangles(struct TVec3<size_t> *tri, size_t N);
void boundary_edges(size_t *boundary, size_t N);
void build_square_mesh(struct MyMesh *m, size_t N, double a = 1.0, double b = 1.0);

void build_disk_vertices(struct TVec2<double> *vert, size_t N, double R = 1.0);
void build_disk_triangles(struct TVec3<size_t> *tri, size_t N);
void disk_boundary_edges(size_t *boundary, size_t N);
void build_disk_mesh(struct MyMesh *m,size_t N, double R = 1.0);


