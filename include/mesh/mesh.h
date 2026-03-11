#pragma once
#include <stdint.h>

#include "array.h"
#include "vec3.h"

#include "hash_table.h"
#include "hash.h"

/* Edge structure */
struct Edge
{
	uint32_t v0;
	uint32_t v1;

	static constexpr uint32_t empty_int = ~uint32_t(0);
	Edge() : v0(empty_int), v1(empty_int) {}
	constexpr Edge(uint32_t a, uint32_t b) : v0(a < b ? a : b), v1(a < b ? b : a) {}

	constexpr bool operator==(const Edge &other) const
	{
		return (v0 == other.v0 && v1 == other.v1);
	}
};

/* Edge structure hasher */
struct EdgeHasher
{
	static constexpr uint32_t empty_int = ~uint32_t(0);
	static constexpr Edge empty_key = Edge(empty_int, empty_int);

	size_t hash(Edge edge) const
	{
		uint32_t hash = 0;
		hash = murmur2_32(hash, edge.v0);
		hash = murmur2_32(hash, edge.v1);
		return hash;
	}

	bool is_empty(Edge edge) const
	{
		return edge == empty_key;
	}

	bool is_equal(Edge edge_1, Edge edge_2) const
	{
		return edge_1 == edge_2;
	}
};

struct Mesh
{
	/* Vertices */
	TArray<Vec3> positions;
	TArray<uint32_t> indices;

	/* Edges */
	TArray<Edge> edges;
	HashTable<Edge, uint32_t, EdgeHasher> edge_idx;

	/* Periodic mappings */
	TArray<uint32_t> periodic_map;
	TArray<uint32_t> dof_map;
	size_t periodic_dofs_count;
	bool is_periodic = false;

	TArray<float> attr;
	size_t vertex_count() const { return positions.size; }
	size_t index_count() const { return indices.size; }
	size_t edge_count() const { return edges.size; }
	size_t triangle_count() const { return indices.size / 3; }

	/* Fills the edge tables in case of P2 FEM */
	void build_edges();
};

inline void Mesh::build_edges()
{
	size_t nt = triangle_count();

	edge_idx.clear();
	edge_idx.reserve(3 * nt);
	uint32_t *current_key;
	uint32_t current_idx = 0;

	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = indices[3 * tri];
		uint32_t b = indices[3 * tri + 1];
		uint32_t c = indices[3 * tri + 2];
		Edge current_edges[3] = {Edge(a, b), Edge(b, c), Edge(c, a)};

		for (size_t k = 0; k < 3; k++)
		{
			current_key = edge_idx.get(current_edges[k]);
			if (!current_key)
			{
				edges.push_back(current_edges[k]);
				edge_idx.set_at(current_edges[k], current_idx++);
			}
		}
	}
}