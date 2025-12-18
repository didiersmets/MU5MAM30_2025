#include <stddef.h>
#include <stdint.h>

#include "array.h"
#include "hash.h"
#include "hash_table.h"
#include "mesh.h"
#include "vec3.h"
#include <iostream>
using namespace std;

/* A proposed solution using a from scratch implementation of a
 * hash table (https://en.wikipedia.org/wiki/Hash_table)
 *
 * The goal of the hash table is to answer
 * the question : given a vertex position (a list of
 * three floating point x,y,z coordinates, here using
 * the type Vec3 defined in common/vec3.h), is that
 * position already in a set of "known" positions ?
 * See include/common/hash_table.h information of use (and implementation)
 *
 * NOTE : Feel free to rely on different implementations like the more
 *        standard std::unordered_map (or even std::map) from the standard
 *        C++ library.
 */

/* A vaudou recipe to turn a (x,y,z) Vec3 into an integer, i.e. a hash map
 * for Vec3. "Good" hash maps should "shuffle" as much as possible the input
 * data to spread it in the table and avoid too many collisions, which recipe
 * is most appropriate is matter of empirism and belief, "murmur_32" is a
 * popular class to hash 4bytes at a time, without any knowledge of what type of
 * data is in these 4 bytes.
 */
struct PositionHasher
{
	const Vec3 *pos;
	static constexpr uint32_t empty_key = ~static_cast<uint32_t>(0);
	size_t hash(uint32_t key) const
	{
		uint32_t hash = 0;
		const uint32_t *p =
			reinterpret_cast<const uint32_t *>(pos + key);
		hash = murmur2_32(hash, p[0]);
		hash = murmur2_32(hash, p[1]);
		hash = murmur2_32(hash, p[2]);
		return hash;
	}

	bool is_empty(uint32_t key) const
	{
		return (key == empty_key);
	}
	bool is_equal(uint32_t key1, uint32_t key2) const
	{
		return pos[key1] == pos[key2];
	}
};

size_t build_position_remap(Vec3 *pos, size_t count, uint32_t *remap)
{
	PositionHasher hasher{pos};
	HashTable<uint32_t, uint32_t, PositionHasher> vtx_remap(count, hasher);

	size_t new_count = 0;
	uint32_t *current_key;

	for (size_t i = 0; i < count; ++i)
	{
		current_key = vtx_remap.get(i);
		if (current_key)
			remap[i] = *current_key;
		else
			remap[i] = new_count++;
		vtx_remap.set_at(i, remap[i]);
	}

	return new_count;
}

void remove_duplicate_vertices(Mesh &m)
{
	Vec3 *pos = m.positions.data;
	size_t vtx_count = m.vertex_count();
	TArray<uint32_t> remap(m.vertex_count());
	size_t index_count = m.index_count();

	/* Build the positions remap */
	size_t new_count = build_position_remap(pos, vtx_count, remap.data);

	TArray<Vec3> new_pos(new_count);
	TArray<uint32_t> new_idx(index_count);

	/* Fill the new positions array */
	for (size_t i = 0; i < vtx_count; i++)
		new_pos[remap[i]] = m.positions[i];

	/* Fill the new indices array */
	for (size_t j = 0; j < index_count; j++)
		new_idx[j] = remap[m.indices[j]];

	m.positions.resize(new_count);
	for (size_t i = 0; i < new_count; ++i)
		pos[i] = new_pos[i];

	for (size_t i = 0; i < index_count; ++i)
		m.indices[i] = new_idx[i];
}
