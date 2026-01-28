#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

static void load_cube_vertices(Vec3 *pos, size_t subdiv);
static void load_cube_indices(uint32_t *idx, size_t subdiv);

int load_cube(Mesh &m, size_t subdiv)
{
	/* Check subdiv is reasonable and return error if not */
	if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */) {
		return (-1);
	}

	size_t n = subdiv + 1;

	/* Reserve memory for vertices and indices */
	m.positions.resize(6 * POW2(n));
	m.indices.resize(36 * POW2(subdiv));

	/* First build vertices as six unattached faces of n^2 vertices each */
	/* See below for implementation */
	load_cube_vertices(m.positions.data, subdiv);

	/* Build corresponding triangulation indices */
	/* See below for implementation */
	load_cube_indices(m.indices.data, subdiv);

	/* Finally attach faces between themselves */
	/* Implementation in src/duplicate_verts.cpp */
	remove_duplicate_vertices(m);

	return (0);
}

static void load_cube_vertices(Vec3 *pos, size_t subdiv)
{
	size_t n = subdiv + 1;
	//size_t h = 1/subdiv;

	///// Commentaires perso/explication - pas vraiment un commentaires de code mais je suis un peu trop habituée à jupyter
	// On va remplir toutes les faces en même temps, mais en faisant attention que les triangles d'une même face se succèdent dans pos (plus facile pour remplir les indices après)
	// On attribu un numéro à chaque face ( front = 0, back = 1, left = 2, right = 3, down = 4, up = 5).
	// on sait que le nombre total de vertices sera 6n**2. donc un multiple de 6. On va faire un remplissage en congruence : quand on remplit
	// un point qui devrait être dans la face k, on accède à pos[6k + i] où i<n**2.
	// pourquoi ça devrait marcher :
			//1) on initialise le tableau pos dans la fonction load cube donc il n'y a pas de problème outre mesure
			//2) on supprime des vertices après avoir construit nos triangles grace à load_indices, donc on peut traiter aussi les faces uen par une pour construire les indices/triangles sans se préocuper de récupérer des sommets avec des numéros bizzares.

	// On va faire une subdivision classique : on a n points par ligne de cadrillage. si on est sur des intervalles [-a;a] on peut prendre le pas
	//h = 2a/subdiv et considérer les points {-a + h*k} avec k = 0,...,subdiv. Quand k = subdiv on a bien -a+h*k = -a+ 2a = a.

	float a = 1.0;
	float h = 2*a/subdiv;
	for(size_t i = 0;i<n; i++){
		float z = -a + h*i;
		for(size_t j = 0;j<n;j++){
			float w = -a + h*j;

			Vec3 v0(w,-a,z); // front // c'est lui qu'on affiche pas je pense !
			pos[j +  i*n] = v0;

			Vec3 v1(w,a,z);
			pos[1*POW2(n) + i*n+ j] = v1; // back


			Vec3 v2(-a,w,z); // left
			pos[2*POW2(n) + i*n+ j]  = v2;

			Vec3 v3(a,w,z); // right
			pos[3*POW2(n) + i*n+ j]  = v3;

			Vec3 v4(w,z,-a); // down
			pos[4*POW2(n) + i*n+ j]  = v4;

			Vec3 v5(w,z,a); // up
			pos[5*POW2(n) +  i*n+j]  = v5;







		}
	}



}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{

	// on reprends les même idées de décallage que précédemment.
	// si on considère i qui parcours une ligne de cadrillage et ligne le numéro de la ligne parcouru. en parcourant de bas en haut et de gauche à droite, le point qui sur le cadrillage se situe en (i,ligne) est numéroté i + lign*n. De cette manière on retrouve les indices des triangles qu'on veut consistuer :
	//
	// 	i + (1 + lign)*n ------ (i+1)+(1+lign)*n
	// 		|						|
	// 		|						|
	// i + lign*n  --------- (i+1) + lign*n


	// Pour retrouve les indices correspondant sur chaque face n°k, il suffit d'ajouter k*n**2. Le remplissage de indices se fait en considérant des modulo 6n**2. En vrai, on va juste faire parcours un compteur et traiter chaque face (avec 6 face donc une petite constante le cout est pas vraiment plus grand) .

	size_t n = subdiv + 1;
	float h = 1/subdiv;

	size_t compteur = 0;

	for(size_t k = 0;k<6;k++){

		for( size_t i = 0; i<n-1; i++){
			for(size_t lign = 0; lign<n-1; lign++){

				// on récupère nos quatre sommets :
				size_t a = i + lign*n ;
				size_t b = (i+1) + lign*n ;
				size_t c = 	i + (1 + lign)*n;
				size_t d = (i+1)+(1+lign)*n;

				// donc on a :
				// c - d
				// | / |
				// a - b



				// on ajoute tout ça à notre liste
				/// premier triangle : acd
				idx[compteur++] = a + k*POW2(n) ;
				idx[compteur++] = c + k*POW2(n) ;
				idx[compteur++] = d+ k*POW2(n) ;
				// deuxième triangle abd
				idx[compteur++] = a + k*POW2(n) ;
				idx[compteur++] = b + k*POW2(n) ;
				idx[compteur++] = d + k*POW2(n) ;
			}
		}
	}








}
