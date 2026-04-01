#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
/* - List of initilisation of references, necessary since
 references cannot be changed, once initilised with a defualt value.
 Necessary, since the constructor is called with a reference &m.
 - N is of type size_t, m.vertex_count() is of type size_t too
 - The rest calls constructors of the respective types (TArray<double>).*/
    : m(m), N(m.vertex_count()), omega(N), Momega(N), psi(N), r(N), p(N), Ap(N)
{
	/* Filling P, M, and S with the right values.*/
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);

	/* Since the sum of the functions phi is 1 on the sphere, the sum
	 of all the matrix entries of the mass matrix is exactly the value
	 of the integral of the 1 fct. over the whole approximated sphere.*/
	vol = M.sum();

	inited = false;
	t = 0;
}


/* Member function */
/* sum(M*V') = sum(M*V) - sum(M*ONES)*s/vol = s-vol*s/vol = s-s = 0 */
void NavierStokesSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i) {
		V[i] -= s / vol; /* V' = V-s/vol*ONES */
	}
}


/* Member function */
void NavierStokesSolver::compute_transport(double *T)
{
	/* T can be an array of type double, since the method is called
	with a pointer anyways. That array is of length N and by applying
	the method memset every byte of that is reserved to staore T will
	be set to 0. */
	memset(T, 0, N * sizeof(double));

	//Rossby number Ro = 0.1 
	double RoRec = 6.28;

	for (size_t t = 0; t < m.triangle_count(); t++) {
		/* indices contains all the indices of points that are contained
		in a trianle.*/
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		assert(a < N && b < N && c < N);
		double sum = omega[a] + omega[b] + omega[c];
		double sumPointsZ = m.positions[a].z + m.positions[b].z + m.positions[c].z;
		/* !!! Calculation needs to be dependent on points on the sphere.
		In this code, it is not. !!! */
		T[a] += (sum + 2*RoRec*sumPointsZ) * (psi[b] - psi[c]);
		T[b] += (sum + 2*RoRec*sumPointsZ) * (psi[c] - psi[a]);
		T[c] += (sum + 2*RoRec*sumPointsZ) * (psi[a] - psi[b]);
	}

	/* The value 1/6 is just the value of \int\phi_{\hat{P}}\,\d x */
	for (size_t v = 0; v < N; v++) {
		T[v] *= 1.0 / 6;
	}
}


/* Member function */
size_t NavierStokesSolver::compute_stream_function()
{
	double b2, r2, rel_error;
	size_t iter;

	/* Defining several pointers, so we do not have to write complicated
	 * references all the time, but can write AP instead of Ap.data, for example.*/
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;
	double *Om = omega.data;
	double *MOm = Momega.data;
	double *Psi = psi.data;

	/* Update M*Omega^m: MOm = M*Om. */
	M.mvp(Om, MOm);

	/* Compute rhs norm2 */
	b2 = blas_dot(MOm, MOm, N);

	/* Form initial R and P */
	S.mvp(Psi, R); // R= S*Psi
	blas_axpby(1, MOm, -1, R, N);
	blas_copy(R, P, N);
	r2 = blas_dot(R, R, N);
	rel_error = sqrt(r2 / b2);

	/*!!! Why do we not use conjugate_gradient_solve? !!!*/

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


/* Member function */
void NavierStokesSolver::time_step(double dt, double nu)
{
	double b2, r2, rel_error;

	size_t iter1, iter2;

	/* Defining several pointers, so we do not have to write complicated
	 * references all the time, but can write AP instead of Ap.data, for example.*/
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
