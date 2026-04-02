#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include "navier_stokes.h"

#include "P1.h"
#include "P2.h"
#include "tiny_blas.h"
#include "cholesky.h"

NavierStokesSolver::NavierStokesSolver(const Mesh &m, int degre, SolverType solver, double nu, double dt)
    : m(m), degre(degre), N(m.vertex_count()), omega(N), Momega(N), psi(N), r(N), p(N), Ap(N), solver_type(solver) {
	if (degre == 1) {
		build_P1_CSRPattern(m, P);
		build_P1_mass_matrix(m, P, M);
		build_P1_stiffness_matrix(m, P, S);
	}
	else if (degre == 2) {
		EdgeAdjacency edge_adj(m);
		build_P2_CSRPattern(m, P, edge_adj);

		N = P.rows; // La nouvelle taille totale
		// On redimensionne tout à la bonne taille car on ne connaissait pas le nombre de DDL avant csr_pattern
		omega.resize(N); Momega.resize(N); psi.resize(N); r.resize(N); p.resize(N); Ap.resize(N);

		// On prépare un tableau de la taille (3 * nombre de triangles)
		edge_ddls.resize(3 * m.triangle_count());

		for (size_t t = 0; t < m.triangle_count(); ++t) {
			uint32_t a = m.indices[3 * t + 0];
			uint32_t b = m.indices[3 * t + 1];
			uint32_t c = m.indices[3 * t + 2];

			// Permet d'avoir un équivalent du mesh mais sur les 3 arêtes qui forment le triangle
			edge_ddls[3 * t + 0] = m.vertex_count() + edge_adj.get_edge_id(a, b);
			edge_ddls[3 * t + 1] = m.vertex_count() + edge_adj.get_edge_id(a, c);
			edge_ddls[3 * t + 2] = m.vertex_count() + edge_adj.get_edge_id(b, c);
		}

		build_P2_mass_matrix(m, P, M, edge_adj);
		build_P2_stiffness_matrix(m, P, S, edge_adj);
	}

	// Initialize Cholesky if requested
	if (solver_type == SolverType::CHOLESKY) {
		/* Create matrix A = M + nu*dt*S */
		// First, store the sparsity pattern
		cholesky_A_pattern.symmetric = S.symmetric;
		cholesky_A_pattern.rows = S.rows;
		cholesky_A_pattern.cols = S.cols;
		cholesky_A_pattern.nnz = S.nnz;
		cholesky_A_pattern.row_start = P.row_start;  // Copy from pattern P
		cholesky_A_pattern.col = P.col;              // Copy from pattern P
		
		//Crée la matrice A avec le pattern stocké et les données calculées
		//Il y avait une erreur dans la gestions des pointeurs qui pointaient vers les données de S et M, 
		//on a du faire une copie pour éviter que les données soient modifiées.
		//Cela causait par ailleurs l'erreur "double free or corruption (out)".
		//Car les données de S et M étaient libérées deux fois : une fois par S et M, et une fois par A.
		//Il faut donc faire une copie des donnees de S et M dans A pour éviter ce problème, et créer une matrice A indépendante de S et M.
		CSRMatrix A;
		A.symmetric = S.symmetric;
		A.nnz = S.nnz;
		A.rows = S.rows;
		A.cols = S.cols;
		A.row_start = cholesky_A_pattern.row_start.data;
		A.col = cholesky_A_pattern.col.data;
		A.data.resize(A.nnz);
		
		// Calcul A.data = M.data + nu*dt*S.data
		for (size_t k = 0; k < A.nnz; k++)
			A.data[k] = M.data[k] + nu * dt * S.data[k];

	
		etree T;
		symbolic_cholesky(A, T);
		L_pattern(A, T, cholesky_L_pattern);
		cholesky_fact(A, cholesky_L, cholesky_L_pattern);
		cholesky_ready = true;
	}

	vol = M.sum();
	inited = false;
	t = 0;
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i) {
		V[i] -= s / vol;
	}
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));
	if (degre == 1) {
		for (size_t t = 0; t < m.triangle_count(); t++) {
			uint32_t a = m.indices[3 * t + 0];
			uint32_t b = m.indices[3 * t + 1];
			uint32_t c = m.indices[3 * t + 2];
			assert(a < N && b < N && c < N);
			double sum = omega[a] + omega[b] + omega[c];
			T[a] += sum * (psi[b] - psi[c]);
			T[b] += sum * (psi[c] - psi[a]);
			T[c] += sum * (psi[a] - psi[b]);
		}
		for (size_t v = 0; v < N; v++) {
			T[v] *= 1.0 / 6;
		}
	}
		else if (degre == 2) { //on va avoir besoin de l'intégration par quadrature de Gauss
struct GaussPoint { double xi, eta, w; };
        const GaussPoint gauss[6] = {
            {0.091576213509771, 0.091576213509771, 0.054975871827661},
            {0.816847572980459, 0.091576213509771, 0.054975871827661},
            {0.091576213509771, 0.816847572980459, 0.054975871827661},
            {0.445948490915965, 0.445948490915965, 0.111690794839005},
            {0.108103018168070, 0.445948490915965, 0.111690794839005},
            {0.445948490915965, 0.108103018168070, 0.111690794839005}
        };

        for (size_t t = 0; t < m.triangle_count(); ++t) {
            uint32_t ddls[6];
            ddls[0] = m.indices[3 * t + 0]; // Sommet A
            ddls[1] = m.indices[3 * t + 1]; // Sommet B
            ddls[2] = m.indices[3 * t + 2]; // Sommet C
            // On a en même temps la correspondance des arêtes grace au tableau
            ddls[3] = edge_ddls[3 * t + 0]; // Arête AB
            ddls[4] = edge_ddls[3 * t + 1]; // Arête AC
            ddls[5] = edge_ddls[3 * t + 2]; // Arête BC
			//On récupère les valeurs (psi et oméga) sur notre triangle
            double psi_loc[6], omega_loc[6];
            for (int i = 0; i < 6; ++i) {
                psi_loc[i]   = psi[ddls[i]];
                omega_loc[i] = omega[ddls[i]];
            }
			//Vecteur local du second membre
            double T_loc[6] = {0.0};
			//ici on intègre sur le triangle car on ne peut plus sommer directement comme en P1
        	//pour cela on boucle sur les 6 points de gauss (car degré 4)
            for (int q = 0; q < 6; ++q) {
                double xi = gauss[q].xi;
                double eta = gauss[q].eta;
                double w = gauss[q].w;

                //Valeur des fonctions de forme (polynômes) en ce point précis (xi, eta)
                double phi[6] = {
                    (1.0 - xi - eta) * (1.0 - 2.0*xi - 2.0*eta),
                    xi * (2.0*xi - 1.0),
                    eta * (2.0*eta - 1.0),
                    4.0 * xi * (1.0 - xi - eta),
                    4.0 * eta * (1.0 - xi - eta),
                    4.0 * xi * eta
                };

                //Dérivées des fonctions de forme
                double dphi_dxi[6] = {
                    -3.0 + 4.0*xi + 4.0*eta,
                     4.0*xi - 1.0,
                     0.0,
                     4.0 - 8.0*xi - 4.0*eta,
                    -4.0*eta,
                     4.0*eta
                };

                double dphi_deta[6] = {
                    -3.0 + 4.0*xi + 4.0*eta,
                     0.0,
                     4.0*eta - 1.0,
                    -4.0*xi,
                     4.0 - 4.0*xi - 8.0*eta,
                     4.0*xi
                };

            	// On utilise les valeurs aux noeuds (psi_loc) et les dérivées (dphi)
            	// pour connaître le gradient exact du fluide à la coordonnée (xi, eta).
                double dpsi_dxi = 0, dpsi_deta = 0, domega_dxi = 0, domega_deta = 0;
                for (int i = 0; i < 6; ++i) {
                    dpsi_dxi    += psi_loc[i]   * dphi_dxi[i];
                    dpsi_deta   += psi_loc[i]   * dphi_deta[i];
                    domega_dxi  += omega_loc[i] * dphi_dxi[i];
                    domega_deta += omega_loc[i] * dphi_deta[i];
                }

                //Jacobien (le transport, terme convectif de navier stokes)
                double transport_val = (dpsi_dxi * domega_deta) - (dpsi_deta * domega_dxi);

            	// L'intégrale est la somme des valeurs calculées, multipliées
            	// par le poids 'w' du point de Gauss et par la fonction test 'phi'.
                for (int i = 0; i < 6; ++i) {
                    T_loc[i] += w * phi[i] * transport_val;
                }
            }

            // Assemblage global classique
            for (int i = 0; i < 6; ++i) {
                T[ddls[i]] += T_loc[i];
            }
        }
		}
}

size_t NavierStokesSolver::compute_stream_function()
{
	double b2, r2, rel_error;
	size_t iter;

	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;
	double *Psi = psi.data;

	M.mvp(Om, MOm);

	/* Compute rhs norm2 */
	b2 = blas_dot(MOm, MOm, N);

	/* Form initial R and P */
	S.mvp(Psi, R);
	blas_axpby(1, MOm, -1, R, N);
	blas_copy(R, P, N);
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	/* Iterate until convergence */
	iter = 0;
	while ((rel_error > tol) && (iter++ < iter_max)) {

		/* Compute AP */
		S.mvp(P, AP);

		/* Update Psi */
		double alpha = r2 / blas_dot(P, AP, N);
		blas_axpy(alpha, P, Psi, N);

		/* Update R */
		blas_axpy(-alpha, AP, R, N);

		/* Update r2 and P */
		double beta = 1.0 / r2;
		r2 = blas_dot(R, R, N);
		rel_error = sqrt(r2 / b2);
		beta *= r2;
		blas_axpby(1, R, beta, P, N);
	}

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	double b2, r2, rel_error;

	size_t iter1, iter2;

	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;

	iter1 = compute_stream_function();

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Form rhs, saved in P */
	compute_transport(P);
	blas_axpby(1, MOm, dt, P, N);

	// Use Cholesky if available
	if (solver_type == SolverType::CHOLESKY && cholesky_ready) {
		TArray<double> omega_perm(N);

		cholesky_solve(cholesky_L, *P, omega_perm);

		M.mvp(Om, MOm);
		set_zero_mean(omega.data);
		t += dt;
		(void)iter1;
	}

	else{
	// Fall back to Conjugate Gradient solver
	b2 = blas_dot(P, P, N);

	/* Form initial R and P */
	S.mvp(Om, R);
	blas_axpby(1, MOm, dt * nu, R, N);
	blas_axpby(1, P, -1, R, N);
	blas_copy(R, P, N);
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	/* Iterate until convergence (and at least once) */
	iter2 = 0;
	do {

		/* Compute AP (invalidates Mom) */
		S.mvp(P, AP);
		M.mvp(P, MOm); /* MOm used as temp storage */
		blas_axpby(1, MOm, dt * nu, AP, N);

		/* Update Om */
		double alpha = r2 / blas_dot(P, AP, N);
		blas_axpy(alpha, P, Om, N);

		/* Update R */
		blas_axpy(-alpha, AP, R, N);

		/* Update r2 and P */
		double beta = 1.0 / r2;
		r2 = blas_dot(R, R, N);
		rel_error = sqrt(r2 / b2);
		beta *= r2;
		blas_axpby(1, R, beta, P, N);

		/* Update MOm */
		M.mvp(Om, MOm);

		iter2++;
	} while ((rel_error > tol) && (iter2 <= iter_max));

	set_zero_mean(omega.data);

	t += dt;

	(void)iter1;
	//printf("Iter 1 : %zu, Iter2 : %zu\n", iter1, iter2);
	}
}