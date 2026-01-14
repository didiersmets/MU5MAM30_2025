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

	size_t new_count = 0; // nombre de sommets uniques

	for (size_t i = 0; i < count; ++i) {
		uint32_t idx;
		if (vtx_remap.find_or_insert(i, &idx)) {
			// sommet déjà présent
			remap[i] = idx;
		} else {
			// sommet nouveau
			remap[i] = new_count;
			vtx_remap.insert(i, new_count);
			new_count++;
		}
	}

	return new_count;
}

void remove_duplicate_vertices(Mesh &m) // Supprime les sommets dupliqués du maillage
{

	Vec3 *pos = m.positions.data; // on récupère le tableau avec les positions doublons
	size_t nb_vtx = m.vertex_count(); //  Nombre initial de sommets

	// On intialise, a priori trop grand, le tableau rempa
	TArray<uint32_t> remap(nb_vtx); // remap[i] = nouvel position du sommet i


	size_t new_count = build_position_remap(pos, nb_vtx, remap.data); // noouveau nombre de sommets, sans duplications


	// changement de positions dans le tableau
	for (size_t i = 0; i < nb_vtx; ++i) { // pour tous les sommest

		// On remplace la position du sommet remap[i] par celle du sommet i
		// on onserve donc la position du sommet unique dans le tableau des positions
		pos[remap[i]] = pos[i];
	}

	// On remet le tableau à la bonne taille
	m.positions.resize(new_count);


	// Changement d'indices '

	// nombre total d'indices dans tous les triangles
	size_t idx_count = m.index_count();

	// Pour chaque indice de chaque triangle
	for (size_t i = 0; i < idx_count; ++i) {
		// On remplace l'ancien index par son index unique
		m.indices[i] = remap[m.indices[i]];
	}

}
