#include "mesh.h"
#include "adjacency.h"

VTAdjacency::VTAdjacency(const Mesh &m)
{
	// On récupère le nombre de noeuds, de triangles et de points à mettre dans vtri = toutes les coordonnées des trois sommets qui composent le triangles
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	size_t idx_count = 3 * tri_count;

	// on alloue la mémoire aux trois tableau de la structure qu'on doit remplir.
	degree.resize(vtx_count);
	offset.resize(vtx_count);
	vtri.resize(idx_count);

	// Il faut initialiser les tableaux degree et offset.
	// Pourquoi pas d'initialisationde vtri ?
	for (size_t i = 0; i < vtx_count; ++i) {
		degree[i] = 0;
		offset[i] = 0;
	}


	for (size_t i = 0; i < idx_count; ++i) {
		// idx_count est la taille de indices  contient tous les sommets, avec répétition.
		// par exemple pour un triangles 0 _ 1
		//       | / |
		// 		 2 _ 3
		// on aurait dans la structure mesh indices = [0,1,2,1,2,3].
		// m.indices[i] accède à chaque élément de indices.
		// et donc si on trouve par exemple i = 2, on ajoute en position d'indice 2 un comptage.
		degree[m.indices[i]]++;
	}

	for (size_t v = 1; v < vtx_count; ++v) {
		offset[v] = offset[v - 1] + degree[v - 1];
	}

	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		vtri[offset[a]++] = { b, c };
		vtri[offset[b]++] = { c, a };
		vtri[offset[c]++] = { a, b };
	}

	for (size_t v = 0; v < vtx_count; ++v) {
		offset[v] -= degree[v];
	}
	assert(offset[vtx_count - 1] + degree[vtx_count - 1] == idx_count);
}
