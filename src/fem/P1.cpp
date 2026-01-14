#include <stdio.h>
#include <string.h>
#include <unordered_set>
#include <algorithm>
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness.h"

/* CSRMatrix variants */

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	// row_start[i] = donne l'indice auquel commence la ligne i dans les tableau col et data. Donc si on a row_start[i] = 3 et row_start[i+1] = 7, on sait que
	// que data[ligne_i] = data[3], data[4], data[5], data[7]'
	// row_start[i+1] - row_start[i] = nombre d’éléments non nuls dans la ligne i
	//col[k] = indice de la colone de la matrice où se trouve le  k-ème élément non nul
	// data[k] = valeur du k-ème coefficient non nul

	size_t nb_vtx = m.vertex_count();
	size_t nb_tri = m.triangle_count();

	//on intialise la strucutre de notre matrice
	P.symmetric = true;
	P.rows = P.cols = nb_vtx;
	P.row_start.resize(nb_vtx + 1);

	// on récupère notre matrice d'adjacence '
	VTAdjacency adj(m);

	size_t max_nnz = 3*nb_tri +nb_vtx; // nombre max de oeffiient non zeros qu'on peut avoir
	// on initialise notre col en fontion '
	P.col.resize(max_nnz);
	// on compte nos non-zeros pour pouvoir ajuster les tailles à la fin (faudrait quand même pas gacher)
	size_t nnz = 0; // nous permet de garder le fil de combien de sommets sont non nul : on sait que quand on a fini une ligne, la prohaine commence excatement à nnz.

	///// dans la suite on va utiliser unordered_set (effacer avant de commit )
	// ens.insert(x)= ajoute x dans le set. Retourne pair<iterator, bool>
	//ens.insert(x).second=bool : true si x a été ajouté, false si x était déjà présent
	//ens.find(x)=renvoie un itérateur vers x si trouvé, ens.end() sinon
	//ens.clear()= vide le set


	for (size_t v = 0; v < nb_vtx; ++v) { // pour chacun des sommets
		P.row_start[v] = nnz; // voir remarque au dessus : on sait qu'on commence notre nouvelle ligne au nnz'
		int32_t start = adj.offset[v]; // on sait que nos données dans adj (donc dans vtri) se trouvent entre offset[i] et offset[i+1] donc entre les indices offset[i] et offset[i] + degree[i]
		uint32_t stop = start + adj.degree[v];

		std::unordered_set<uint32_t> seen; // on crée notre ensemble non ordonné d'éléments déjà vus


		for (size_t k = start; k < stop; ++k) { // on considère le triangle u-v-w
			uint32_t w = adj.vtri[k].next;
			uint32_t u = adj.vtri[k].prev;

			bool not_seen_w = seen.insert(w).second; // est ce qu'on a déjà insérer w ? '=> insert essaie d'ajouter x. Faire insert.second c'est répondre à la question : si j'essaie de l'ajouter, est ce que ça fontionne ? si ça ne fonctionne pas, c'est qu'on l'avait déjà dedans => false = x était déjà là.
			bool not_seen_u = seen.insert(u).second; // est ce qu'on a déjà insérer u ? '

			if (w < v && not_seen_w) {
				P.col[nnz++] = w;
				seen.insert(w);
			}

			if (u < v && not_seen_u) {
				P.col[nnz++] = u;
				seen.insert(u);
		}
		P.col[nnz++] = v;
	}}

	P.row_start[nb_vtx] = nnz;
	P.col.resize(nnz);
	P.col.shrink_to_fit();

	for (size_t i = 0; i < nb_vtx; ++i) {

		uint32_t* line_start = &P.col[P.row_start[i]];// pointeur vers le début de la ligne i dans P.col
		size_t count = P.row_start[i + 1] - P.row_start[i]; // nombre d'éléments dans cette ligne

		std::sort(line_start, line_start + count);// tri croissant des indices de colonnes de cette ligne

	}
}
void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	assert(P.row_start.size == vtx_count + 1);

	M.symmetric = true;
	M.rows = M.cols = vtx_count;
	M.nnz = P.col.size;
	M.row_start = P.row_start.data;
	M.col = P.col.data;
	M.data.resize(M.nnz);
	for (size_t i = 0; i < M.nnz; ++i) {
		M.data[i] = 0.0;
	}

	// petite fonction lambda utile pour la suite :
	// on ajoute que en faisan bien attention à l'ordre des idnices ;'
	auto add_offdiag = [&](uint32_t i, uint32_t j, double val){
		if (i < j) std::swap(i, j); // i >= j
		M(i,j) += val;
	};


	const TArray<uint32_t> &idx = m.indices;

	for (size_t i = 0; i < tri_count; ++i) {
		// on récupère les sommets et leurs coordonnées
		uint32_t a = m.indices[3*i]; Vec3f A = m.positions[a];
		uint32_t b = m.indices[3*i+1]; Vec3f B = m.positions[b];
		uint32_t c = m.indices[3*i+2]; Vec3f C = m.positions[c];

		// on construit les vecteurs AB et AC qu'on utilise pour calculer l
		Vec3d AB = Vec3d(B) - Vec3d(A);
		Vec3d AC = Vec3d(C) - Vec3d(A);

		// on initialise notre matrice de masse locale. Si j'ai bien compris on ne doit pas utiliser le format de matrices fem, donc on prends juste un tableau de deux valeurs (0 = coef diag, 1 = coef off diag )'
		double Mloc[2];

		// et on calcul la matrice de masse locale grâce à la fonction de mass.h
		mass(AB, AC, Mloc);

		// maintenant, il faut remplir la matrice de masse globale. (voir dessin cahier mais en gros on ajoute notre notre coef local (i,j) aau coef qui correspond dans la matrice globale)
		// pour les diagonales, pas de problèmes
		M(a, a) += Mloc[0];
		M(b, b) += Mloc[0];
		M(c, c) += Mloc[0];
		// puis dans la boucle triangle
		add_offdiag(a, b, Mloc[1]);
		add_offdiag(b, c, Mloc[1]);
		add_offdiag(c, a, Mloc[1]);

	}
}
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();
	assert(P.row_start.size == vtx_count + 1);

	// initialisation CSR
	S.symmetric = true;
	S.rows = S.cols = vtx_count;
	S.nnz = P.col.size;
	S.row_start = P.row_start.data;
	S.col = P.col.data;
	S.data.resize(S.nnz);
	for (size_t i = 0; i < S.nnz; ++i) S.data[i] = 0.0;

	const TArray<uint32_t> &idx = m.indices;

	// lambda pour ajouter off-diagonale en respectant la symétrie
	auto add_offdiag = [&](uint32_t i, uint32_t j, double val){
		if (i < j) std::swap(i, j); // i >= j
		S(i,j) += val;
	};


	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = idx[3*t + 0];
		uint32_t b = idx[3*t + 1];
		uint32_t c = idx[3*t + 2];

		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];

		Vec3d AB = Vec3d(B) - Vec3d(A);
		Vec3d AC = Vec3d(C) - Vec3d(A);

		// calcul de la matrice locale de rigidité
		double Sloc[6];
		stiffness(AB, AC, Sloc);

		// diagonale
		S(a,a) += Sloc[0];
		S(b,b) += Sloc[1];
		S(c,c) += Sloc[2];

		// off-diagonale
		add_offdiag(a,b, Sloc[3]);
		add_offdiag(b,c, Sloc[4]);
		add_offdiag(c,a, Sloc[5]);
	}
}


/* FEMatrix variants */

void build_P1_mass_matrix(const Mesh &m, FEMatrix &M)
{
	size_t vtx_count = m.vertex_count();
	size_t tri_count = m.triangle_count();

	M.fem_type = FEMatrix::P1_cst;
	M.m = &m;
	M.rows = M.cols = vtx_count;

	M.diag.resize(vtx_count);
	memset(M.diag.data, 0, vtx_count * sizeof(double));

	M.off_diag.resize(tri_count);
	const TArray<uint32_t> &idx = m.indices;
	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = idx[3 * t + 0];
		uint32_t b = idx[3 * t + 1];
		uint32_t c = idx[3 * t + 2];
		Vec3f A = m.positions[a];
		Vec3f B = m.positions[b];
		Vec3f C = m.positions[c];
		Vec3d AB = { (double)B[0] - (double)A[0],
			     (double)B[1] - (double)A[1],
			     (double)B[2] - (double)A[2] };
		Vec3d AC = { (double)C[0] - (double)A[0],
			     (double)C[1] - (double)A[1],
			     (double)C[2] - (double)A[2] };
		double Mloc[2];
		mass(AB, AC, Mloc);
		M.diag[a] += Mloc[0];
		M.diag[b] += Mloc[0];
		M.diag[c] += Mloc[0];
		M.off_diag[t] = Mloc[1];
	}
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();
    assert(P.row_start.size == vtx_count + 1);

    // initialisation CSR
    S.symmetric = true;
    S.rows = S.cols = vtx_count;
    S.nnz = P.col.size;
    S.row_start = P.row_start.data;
    S.col = P.col.data;
    S.data.resize(S.nnz);
    for (size_t i = 0; i < S.nnz; ++i) S.data[i] = 0.0;

    const TArray<uint32_t> &idx = m.indices;

    // lambda pour ajouter off-diagonale en respectant la symétrie
    auto add_offdiag = [&](uint32_t i, uint32_t j, double val){
        if (i < j) std::swap(i,j);
        S(i,j) += val;
    };

    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t + 0];
        uint32_t b = idx[3*t + 1];
        uint32_t c = idx[3*t + 2];

        Vec3f A = m.positions[a];
        Vec3f B = m.positions[b];
        Vec3f C = m.positions[c];

        Vec3d AB = Vec3d(B) - Vec3d(A);
        Vec3d AC = Vec3d(C) - Vec3d(A);

        // calcul de la matrice locale de rigidité
        double Sloc[6];
        stiffness(AB, AC, Sloc);

        // diagonale
        S(a,a) += Sloc[0];
        S(b,b) += Sloc[1];
        S(c,c) += Sloc[2];

        // off-diagonale
        add_offdiag(a,b, Sloc[3]);
        add_offdiag(b,c, Sloc[4]);
        add_offdiag(c,a, Sloc[5]);
    }
}
