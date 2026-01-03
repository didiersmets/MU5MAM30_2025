#include <stddef.h>
#include <stdint.h>

#include "array.h"
#include "hash.h"
#include "hash_table.h"
#include "mesh.h"
#include "vec3.h"

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
struct PositionHasher {
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
	PositionHasher hasher{ pos };
	HashTable<uint32_t, uint32_t, PositionHasher> vtx_remap(count, hasher);

	/* Your implementation goes here */

	/*
	what we want to do is for each vertex to check if the same vertex is already
	in the hash table:
	if it is we set remap[i] to the index of the already existing vertex
	if it is not we add it to the hash table and set remap[i] to the current new index
	and we increment the new vertex count to keep track of the number of unique vertices 
	that have been added so far
	*/

	size_t new_count = 0;

	for( int i = 0; i < count; i++){
		//find if the vertex is already in the hash table
		uint32_t *vtx = vtx_remap.get_or_set(i, new_count);
		//if it is not in the hash table, we add it and set remap[i] to new_count
		remap[i] = (vtx) ? *vtx : new_count++;
	}

	return new_count;
}

void remove_duplicate_vertices(Mesh &m)
{
	Vec3 *pos = m.positions.data;
	size_t vtx_count = m.vertex_count();
	TArray<uint32_t> remap(m.vertex_count());

	size_t new_count = build_position_remap(pos, vtx_count, remap.data);

	/* Remap vertices */
	/* Your implementation goes here */

	/*
	now that we know the number of added vertices and the new remap
	we can actually remove the duplicate vertices by overwriting the
	positions array with the unique vertices only
	the main idea is that:

	remap[i] tells us the new index for vertex i

	If vertex i is the first occurrence of its position --> remap[i] == i
	for unique vertices we writes to the same place 

	If vertex i is a duplicate --> remap[i] < i, pointing to the index of the first occurrence
	for duplicates we overwrite the first occurrence with the same value 

	*/

	for( int i = 0; i < vtx_count; i++)
	 	pos[remap[i]] = pos[i];
	m.positions.resize(new_count);

	/* Remap indices */
	/* Your implementation goes here */

	/*
	we successfully removed the duplicate vertices but the triangles still reference
	the old vertex indices, so we need to remap them using the remap array
	the idea is to map each triangle index as its remap value in its own position

	triangle_index[i] = remap[triangle_index[i]]

	*/

	for (int i = 0; i < m.index_count(); i++)
    	m.indices[i] = remap[m.indices[i]];

}
