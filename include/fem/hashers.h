#include "hash_table.h"
#include "hash.h"
#include "mesh.h"

/* Vertex pairs structure */
struct VertexPair
{
    uint32_t first;
    uint32_t second;
    bool operator==(const VertexPair &other) const
    {
        return first == other.first && second == other.second;
    }
};

/* Vertex pairs structure hasher */
struct VertexPairHasher
{
    static constexpr uint32_t empty_int = ~uint32_t(0);
    static constexpr VertexPair empty_key{empty_int, empty_int};

    size_t hash(VertexPair vertex_pair) const
    {
        uint32_t hash = 0;
        hash = murmur2_32(hash, vertex_pair.first);
        hash = murmur2_32(hash, vertex_pair.second);
        return hash;
    }

    bool is_empty(VertexPair vertex_pair) const
    {
        return vertex_pair.first == empty_int && vertex_pair.second == empty_int;
    }

    bool is_equal(VertexPair vertex_pair_1, VertexPair vertex_pair_2) const
    {
        return vertex_pair_1.first == vertex_pair_2.first && vertex_pair_1.second == vertex_pair_2.second;
    }
};

/* Edge pairs structure */
struct EdgePair
{
    Edge first;
    Edge second;

    constexpr EdgePair() : first(0, 0), second(0, 0) {}
    constexpr EdgePair(const Edge &edge_1, const Edge &edge_2) : first(edge_1), second(edge_2) {}

    constexpr bool operator==(const EdgePair &other) const
    {
        return first == other.first && second == other.second;
    }
};

/* Edge pairs structure hasher */
struct EdgePairHasher
{
    static constexpr uint32_t empty_int = ~uint32_t(0);
    static constexpr EdgePair empty_key{Edge(empty_int, empty_int), Edge(empty_int, empty_int)};

    size_t hash(const EdgePair &edge_pair) const
    {
        uint32_t hash = 0;
        hash = murmur2_32(hash, edge_pair.first.v0);
        hash = murmur2_32(hash, edge_pair.first.v1);
        hash = murmur2_32(hash, edge_pair.second.v0);
        hash = murmur2_32(hash, edge_pair.second.v1);
        return hash;
    }

    bool is_empty(const EdgePair &edge_pair) const
    {
        return edge_pair == empty_key;
    }

    bool is_equal(const EdgePair &edge_pair_1, const EdgePair &edge_pair_2) const
    {
        return edge_pair_1 == edge_pair_2;
    }
};

/* Vertex-Edge pairs structure */
struct VertexEdgePair
{
    Edge edge;
    uint32_t vertex_idx;

    constexpr VertexEdgePair() : edge(0, 0), vertex_idx(0) {}
    constexpr VertexEdgePair(const Edge &e, const uint32_t &idx) : edge(e), vertex_idx(idx) {}

    constexpr bool operator==(const VertexEdgePair &other) const
    {
        return edge == other.edge && vertex_idx == other.vertex_idx;
    }
};

/* Vertex-Edge pairs structure hasher */
struct VertexEdgePairHasher
{
    static constexpr uint32_t empty_int = ~uint32_t(0);
    static constexpr VertexEdgePair empty_key{Edge(empty_int, empty_int), empty_int};

    size_t hash(const VertexEdgePair &vertex_edge_pair) const
    {
        uint32_t hash = 0;
        hash = murmur2_32(hash, vertex_edge_pair.edge.v0);
        hash = murmur2_32(hash, vertex_edge_pair.edge.v1);
        hash = murmur2_32(hash, vertex_edge_pair.vertex_idx);
        return hash;
    }

    bool is_empty(const VertexEdgePair &vertex_edge_pair) const
    {
        return vertex_edge_pair == empty_key;
    }

    bool is_equal(const VertexEdgePair &vertex_edge_pair_1, const VertexEdgePair &vertex_edge_pair_2) const
    {
        return vertex_edge_pair_1 == vertex_edge_pair_2;
    }
};
