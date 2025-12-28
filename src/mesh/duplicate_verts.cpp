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

	
	size_t new_count=0;
	for (size_t i=0;i<count;i++){
		/*i check if pos[i] is already mapped in vtx_map, */
		if (vtx_remap.get(pos[i])==nullptr){ /*if not present*/
			vtx_remap.set(i, new_count); /*map pos[i] to new_count*/
			remap[i]=new_count; /*remap associa al suo i-esimo elemnto il valore che assume nella hash map*/
			new_count++; /*increment new_count since I have added a new unique vertex*/
		}
		else{
			remap[i]=vtx_remap.get(pos[i]); 	/*remap[i] is the index of the already existing vertex*/
		}
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
	m.positions.resize(new_count);
	std::vector<bool> copied(new_count, false);/*I create a vector which tracks which positions have been copied, initialized to false*/
	/* Copy the unique verticeses*/
	for (size_t i=0;i<vtx_count;i++){
		if (!copied[remap[i]]){ /*if the position has not been copied yet*/		
			m.positions[remap[i]]=pos[i];
			copied[remap[i]]=true;
		}
	}
	
	/* Remap indices */
	for (size_t i = 0; i < m.indices.size(); ++i) {
        m.indices[i] = remap[m.indices[i]];
    }
	/* Your implementation goes here */
}
