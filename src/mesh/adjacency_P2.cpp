//
// Created by aurel on 05/03/2026.
//

#include "../../include/mesh/adjacency_P2.h"

EdgeAdjacency::EdgeAdjacency(const Mesh &m) {
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();
    size_t idx_count = 3 * tri_count;
    degree.resize(vtx_count);
    offset.resize(vtx_count);

    for (size_t i = 0; i < vtx_count; ++i) {
        degree[i] = 0;
        offset[i] = 0;
    }
    /*on ajoute 2 degrés cette fois-ci car on stocke b et c séparement*/
    for (size_t i = 0; i < idx_count; ++i) {
        degree[m.indices[i]] = degree[m.indices[i]] + 2;
    }

    for (size_t v = 1; v < vtx_count; ++v) {
        offset[v] = offset[v - 1] + degree[v - 1];
    }

    // Allocation du tableau des arêtes (avec les doublons pour l'instant)
    edges.resize(offset[vtx_count - 1] + degree[vtx_count - 1]);

    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = m.indices[3 * t + 0];
        uint32_t b = m.indices[3 * t + 1];
        uint32_t c = m.indices[3 * t + 2];
        edges[offset[a]++] = { b, 0 };
        edges[offset[a]++] = { c, 0 };
        edges[offset[b]++] = { c, 0 };
        edges[offset[b]++] = { a, 0 };
        edges[offset[c]++] = { a, 0 };
        edges[offset[c]++] = { b, 0 };
    }
/*on note tout les voisins à la manière d'adjacency*/
    for (size_t v = 0; v < vtx_count; ++v) {
        offset[v] -= degree[v];
    } /*on remet les offset comme dans adjacency

    A cette étape on a un grand tableau avec tout les voisins mais plein de doublons
    dans P1 on n'en tenait pas compte et on s'en occupait avec la fonction find pendant la création
    de CSR_pattern mais c'était possible car les sommets étaient déjà indicés de manière unique.
    Ici on n'a 2 problèmes, on ne connait pas le nombre d'arêtes donc on ne va pas avoir la taille de CSR
    à l'avance et on n'a pas de numéro unique d'arête ce qui va poser soucis également pour ajouter l'arête
    dans la matrice lorsque 2 triangles en ont une en commun. On décide donc de faire la table d'adjacence complète
    dans une fonction à part après la création de cette matrice dans P2.cpp

    On va donc faire le tri ici comme pour P1 puis les doublons seront collés donc on va les supprimer puis
    mettre à jour degré*/

    for (size_t a = 0; a < vtx_count; ++a) {
        /* "__restrict" technique d'optimisation venant de P1*/
        EdgeRecord *__restrict to_sort = &edges[offset[a]]; /* "&" pointeur vers l'adresse car sinon fait une copie*/
        /* Tri par insertion comme dans P1 */
        for (size_t k = 1; k < degree[a]; ++k) {
            size_t j = k;
            while (j && (to_sort[j - 1].neighbor > to_sort[j].neighbor)) {
                EdgeRecord tmp = to_sort[j - 1];
                to_sort[j - 1] = to_sort[j];
                to_sort[j] = tmp;
                j--;
            }
        }

        size_t write_idx = 1;

        for (size_t read_idx = 1; read_idx < degree[a]; ++read_idx) {

            if (to_sort[read_idx].neighbor != to_sort[write_idx - 1].neighbor) {
                to_sort[write_idx] = to_sort[read_idx];
                write_idx++;
            }
        }
        degree[a] = write_idx;
    }
    /*maintenant on a plus de doublons mais les offset sont faux et le tableau est trop grand, il faut
     *également numéroter les arêtes, on va faire tout cela d'un coup*/

    size_t write_pos = 0;
    num_edges = 0;
    for (size_t a = 0; a < vtx_count; ++a) {
        size_t old_offset = offset[a]; /*on garde l'ancien pour naviguer dans le tableau*/
        offset[a] = write_pos; /*on mets à jour offset*/
        for (uint32_t k = 0; k < degree[a]; ++k) {
            uint32_t b = edges[old_offset + k].neighbor;
            if (a < b) { /*on n'a pas encore croisé l'arête*/
                edges[write_pos++] = { b, num_edges };
                num_edges++;
            }
            /*sinon on l'a déjà vu donc le degré et le offset est bon car il a déjà été changé,
             *on assigne juste l'identifiant de l'arête)*/
            else if (a > b) {
                edges[write_pos++] = { b, get_edge_id(b, a) };
            }
        }
    }
    edges.resize(write_pos);
}
//Fonction pratique qui pour 2 sommets renvoie l'identifiant de l'arête
uint32_t EdgeAdjacency::get_edge_id(uint32_t a, uint32_t b) const {
    // On récupère le point de départ et le nombre de voisins du sommet 'a'
    size_t start = offset[a];
    size_t count = degree[a];

    // On parcourt uniquement la petite liste des voisins de 'a'
    for (size_t i = 0; i < count; ++i) {
        if (edges[start + i].neighbor == b) {
            return edges[start + i].edge_id; //on l'a trouvé donc on renvoie l'identifiant
        }
    }
    return 0;
}