#pragma once

#include <cstdio>
#include "mesh.h"
#include <vector>

// Print a std::vector of generic type
template <typename T>
inline void print_vector(const std::vector<T>& vec, const char* name = "Vector") {
	printf("%s: [", name);
	for (size_t i = 0; i < vec.size(); ++i) {
		printf("%lu", vec[i]);
		if (i < vec.size() - 1) {
			printf(", ");
		}
	}
	printf("]\n");
}

// Print a std::vector as triplets (for indices representing triangles)
template <typename T>
inline void print_vector_triplets(const std::vector<T>& vec, const char* name = "Indices") {
	printf("%s (as triangles): ", name);
	for (size_t i = 0; i < vec.size(); i += 3) {
		printf("(%lu, %lu, %lu)", vec[i], vec[i + 1], vec[i + 2]);
		if (i + 3 < vec.size()) {
			printf("  ");
		}
	}
	printf("\n");
}

template <typename T>
void print_point(const TVec3<T>& p) {
	printf("Point: (x,y,z) = (%lf, %lf, %lf)\n", p.x, p.y, p.z);
}

// Print mesh information: number of vertices, triangles, etc.
inline void print_mesh_info(const Mesh& mesh) {
	printf("=== Mesh Information ===\n");
	printf("Vertex count: %zu\n", mesh.vertex_count());
	printf("Triangle count: %zu\n", mesh.triangle_count());
	printf("Index count: %zu\n", mesh.index_count());
}

// Print a single face of the mesh
// vertex_start: starting index for vertices to print
// vertex_end: ending index for vertices to print (exclusive)
// index_start: starting index for triangles to print
// index_end: ending index for triangles to print (exclusive)
// row_size: number of vertices per row
inline void print_face(const Mesh& mesh, size_t vertex_start, size_t vertex_end, 
                      size_t index_start, size_t index_end, size_t row_size = 5) {
	printf("\n=== Face Vertices ===\n");
	for (size_t i = vertex_start; i < vertex_end; ++i) {
		printf("(%f, %f, %f)", mesh.positions.data[i].x, mesh.positions.data[i].y, mesh.positions.data[i].z);
		
		if ((i - vertex_start + 1) % row_size == 0) {
			printf("\n");
		} else if (i < vertex_end - 1) {
			printf("  ");
		}
	}
	if ((vertex_end - vertex_start) % row_size != 0) {
		printf("\n");
	}
	
	printf("\n=== Face Indices (as triangles) ===\n");
	for (size_t i = index_start; i < index_end; i += 3) {
		printf("(%u, %u, %u)", mesh.indices.data[i], mesh.indices.data[i + 1], mesh.indices.data[i + 2]);
		if (i + 3 < index_end) {
			printf("  ");
		}
	}
	printf("\n");
}

// Print cube mesh faces
// num_faces: number of faces to print (1-6). Faces are printed in order: UP, DOWN, FRONT, RIGHT, REAR, LEFT
// row_size: number of vertices per row
inline void print_square_mesh(const Mesh& mesh, size_t num_faces = 6, size_t row_size = 5) {
	// Validate num_faces
	if (num_faces == 0 || num_faces > 6) {
		printf("Error: num_faces must be between 1 and 6\n");
		return;
	}
	
	size_t vertices_per_face = mesh.vertex_count() / num_faces;
	size_t indices_per_face = mesh.index_count() / num_faces;
	
	printf("\n=== Square Mesh (%zu Faces) ===\n", num_faces);
	printf("Total vertices: %zu | Per face: %zu | mesh.vertex_count() \n", mesh.vertex_count(), vertices_per_face);
	printf("Total indices: %zu | Per face: %zu\n\n", mesh.index_count(), indices_per_face);
	
	const char* face_names[] = {"UP", "DOWN", "FRONT", "RIGHT", "REAR", "LEFT"};
	
	for (size_t face = 0; face < num_faces; ++face) {
		printf("===== Face %zu: %s =====\n", face, face_names[face]);
		
		size_t v_start = face * vertices_per_face;
		size_t v_end = (face + 1) * vertices_per_face;
		size_t i_start = face * indices_per_face;
		size_t i_end = (face + 1) * indices_per_face;
		
		print_face(mesh, v_start, v_end, i_start, i_end, row_size);
	}
}
