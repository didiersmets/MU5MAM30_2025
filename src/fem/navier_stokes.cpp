#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

//
#include "adjacency.h"
#include "mesh.h"
#include "conjugate_gradient.h"
#include "vec3.h"



NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
	, r(N)
	, p(N)
	, Ap(N)
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);
#endif
	vol = M.sum();
	inited = false;
	t = 0;
}

void NavierStokesSolver::set_zero_mean(double *V)
{

		double mean = 0.0;

		for (size_t i = 0; i < N; ++i)
			mean += M(i, i)*V[i];

		mean /= vol;

		for (size_t i = 0; i < N; ++i)
			V[i] -= mean;


}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	// T(\Omega, \Psi) = \sum_{triangles} \sum_{i,j } \Omega_i \Psi_j (\int_triangle \phi_i \Nabla \phi_j ^T \nabla \phi_k )

	size_t nb_tri = m.triangle_count();
	// on récupère notre matrice d'adjacence '
	VTAdjacency adj(m); // la liste des triangles = adj.vtri[k]


	for(size_t t = 0; i<nb_tri ; i++ ){

		size_t areaT = norm(cross((s0 - s1),(s1-s2)))/2  // |(s_0 - s1) x (s_1 - s_2)|/2

		//les indies
		size_t i0 = mesh.indices[3*t];
		size_t i1 = mesh.indices[3*t + 1];
		size_t i2 = mesh.indices[3*t + 2 ];

		// les sommets en coordonnées :
		vec3 s0 = mesh.positions[i0];
		vec3 s1 = mesh.positions[i1];
		vec3 s2 = mesh.positions[i2];
		TArray<vec3> sommets(s0,s1,s2);

		// les normales : // formules du cours de M1
		vec3 n0 = cross((s1-s2),(s1-s0));
		vec3 n1 = cross((s2-s0),(s2-s1));
		vec3 n2 = cross((s0-s1),(s0-s2));

		TArray<vec3> normales(n0,n1,n2);
		// les gradient des fontions de formes


		vec3 grad_i0 = n0 / dot((s0-s1),n0);
		vec3 grad_i1 = n1 / dot((s1-s2),n1);
		vec3 grad_i2 = n2 / dot((s2-s0),n2);

		TArray<vec3> grad(grad_i0,grad_i1,grad_i2);

		TArray<double>	Om(omega[i0],Omega[i1],Omega[i2]);
		TArray<double> Psi(psi[i0],psi[i1],psi[i2]);
		// on pourrait faire une boucle pour les trois cas (i,j,k) mais je sais pas comment faire et ça ira plus vite de taper que de trouver




		for(size_t i = 0; i<3; i++){
			for(size_t j = 0; j<3; j++){

			double prod = Om[i]*Psi[j];

				for(size_t k = 0; k<3; k++){


					T += prod*areaT/3*dot(grad[j],grad[k])

				}
			}
		}

	}

}

size_t NavierStokesSolver::compute_stream_function()
{
	// initialisation ?
	size_t iter = 0;
	//size_t max_iter = 10*N; il est déjà déclaré dans navier_stokes.h

	//memset(psi.data, 0, N * sizeof(double)); Non parce qu'on s'en sert mis à jours

	// On fait une résolution complète
	memset(r.data, 0, N * sizeof(double));
	memset(p.data, 0, N * sizeof(double));
	memset(Ap.data, 0, N * sizeof(double));





	double rel_error, b2, ;
	b2 = blas_dot(Momega, Momega, N);
	r2 = blas_dot(r,r, N);
	rel_error = sqrt(r2 / b2);


	// résolution par gradient conjugué

	if (!inited) {
		/* r_0 = b - Ax_0 */
		S.mvp(psi, r);
		blas_axpby(1, Momega, -1, r, N);
		/* p_0 = r_0 */
		blas_copy(r, p, N);
	}
	while ((iter < max_iter) && (*rel_error > tol))
	{
		const size_t N = S.rows;
		S.mvp(p, Ap);

		//alpha k-1
		double pAp = blas_dot(p, Ap, N);//<p_{k-1}, Ap_{k-1}>
		double alpha = r2 / pAp;

		// mise à jour de x => x<- x + alpha*p (on a une forme axpy )
		blas_axpy(alpha, p, , N);

		// mise à jour de r => r<- r - alpha Ap
		blas_axpy(-alpha, Ap, r, N);

		//on calcul le nouveau r², on a toujours le ||r_{k-1}||**2 en stock dans r2
		double r2_new = blas_dot(r, r, N);

		// on calcul beta qui est  ||r_k||**2 / ||r_{k-1}||**2
		double beta = r2_new / r2;

		// On met à jorus pk => p<-r + beta*p
		blas_axpby(1.0, r, beta, p, N);

		iter++;

	}
	//(S, // A
										  // Momega,//b
										//psi, //x, Ax=b

	// le fameux problème d'infinité de solution avec psi+C'

	set_zero_mean(psi.data);




	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	compute_stream_function(); // On compute tout par référence S*Psi^{m+1} = -M*\Omega^m. On le fait en premier par ce quand on va
								// résoudre (M+nu*dt*S)*Omega = b (avec b qui fait intervenir Omega^m) on obtient Omega^{m+1}

	// On a besoin de T(Omega^m,psi^m) => Psi est mit à jours par référence donc on sait que c'est bien psi^m avec lequel on travaille '
	//TArray<double> T(N); => on a pas besoin d'initialiser parce que dans compute_transport on met T à 0
	compute_transport(*T);

	TArray<vec3> b = Momega + dt*T ;

	// Maintenant il faut résoudre b = omega*(M + nu*dt*S)


	// On fait une résolution complète donc on initialise
	memset(r.data, 0, N * sizeof(double));
	memset(p.data, 0, N * sizeof(double));
	memset(Ap.data, 0, N * sizeof(double));



	//size_t omega = conjugate_gradient_solve(M + nu*dt*S, // A
										  // b,//b
										//omega, //x, Ax=b
									//	r, // résidus r = Ax-b
									//	p, // initi
									//	Ap,// initi
									//	&rel_error,// à init
									//	tol,
									//	iter_max, // initi
										//inited);

	CSRMatrix A = M + nu*dt*S;
	double b2 = blas_dot(b, b, N);

	if (!inited) {
		/* r_0 = b - Ax_0 */
		A.mvp(x, r);
		blas_axpby(1, b, -1, r, N);
		/* p_0 = r_0 */
		blas_copy(r, p, N);
	}

	double r2 = blas_dot(r, r, N);
	*rel_error = sqrt(r2 / b2);

	int iter = 0;
	while ((iter < max_iter) && (*rel_error > tol)) {
		const size_t N = A.rows;
		A.mvp(p, Ap);

		//alpha k-1
		double pAp = blas_dot(p, Ap, N);//<p_{k-1}, Ap_{k-1}>
		double alpha = r2 / pAp;

		// mise à jour de x => x<- x + alpha*p (on a une forme axpy )
		blas_axpy(alpha, p, x, N);

		// mise à jour de r => r<- r - alpha Ap
		blas_axpy(-alpha, Ap, r, N);

		//on calcul le nouveau r², on a toujours le ||r_{k-1}||**2 en stock dans r2
		double r2_new = blas_dot(r, r, N);

		// on calcul beta qui est  ||r_k||**2 / ||r_{k-1}||**2
		double beta = r2_new / r2;

		// On met à jorus pk => p<-r + beta*p
		blas_axpby(1.0, r, beta, p, N);

		*rel_error = sqrt(r2 / b2);

		iter++;
	}



	// on met à jours Momega :
	Momega = M*omega;

	set_zero_mean(omega.data);

	t += dt;
}
