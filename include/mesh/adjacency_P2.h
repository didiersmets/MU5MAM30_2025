//
// Created by aurel on 05/03/2026.
//

#ifndef MU5MAM30_ADJACENCY_P2_H
#define MU5MAM30_ADJACENCY_P2_H
#include <stdint.h>
#include "array.h"
#include "mesh.h"

/*****************************************************************************
 * Pour P2, on va chercher à indicer chaque arête du maillage de manière unique
 * comme pour les sommets.
 * * EdgeAdjacency construit un graphe d'adjacence Sommet-à-Sommet sans doublon,
 * et attribue un identifiant unique (edge_id) à chaque arête.
 *
 * Une fois chargée, la structure contient :
 *
 * - num_edges : le nombre total d'arêtes uniques dans le maillage.
 * - degree    : un tableau de taille m.vertex_count(). degree[a] est le nombre
 * de voisins (arêtes connectées) du sommet a.
 * - offset    : un tableau de taille m.vertex_count(). Indique la position de
 * départ dans le tableau 'edges' pour le sommet a on y navigue comme pour adjacency :
 * offset[0] = 0 et offset[k] = offset[k-1] + degree[k-1]
 * - edges     : un grand tableau plat qui contient les voisins et l'identifiant de l'arête.
 * Les voisins du sommet k se trouvent aux indices j tels que
 * offset[k] <= j < offset[k] + degree[k].
 *
 *****************************************************************************/

struct EdgeAdjacency {
    struct EdgeRecord {
        uint32_t neighbor; // Le sommet voisin
        uint32_t edge_id;  // L'identifiant unique de l'arête entre le sommet et 'neighbor'
    };

    uint32_t num_edges;
    TArray<uint32_t> degree;
    TArray<uint32_t> offset;
    TArray<EdgeRecord> edges;

    /* Le constructeur qui est implémenté dans src/adjacency_P2.cpp */
    EdgeAdjacency(const Mesh &m);

    /* * Fonction utilitaire pour récupérer l'identifiant de l'arête (a, b).
     */
    uint32_t get_edge_id(uint32_t a, uint32_t b) const;
};

#endif
