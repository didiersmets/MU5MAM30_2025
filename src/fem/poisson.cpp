#include "poisson.h"
#include "square.h"
#include "P1.h"
#include "array.h"
#include "conjugate_gradient.h"
#if USE_FEM_MATRIX
#include "fem_matrix.h"
#else
	#include "sparse_matrix.h"
#endif
#include "mesh.h"
#include "tiny_blas.h"
#include "boundary.h"

PoissonSolver::PoissonSolver(const Mesh &m)
    : m(m), N(m.vertex_count()), f(N), u(N, 0.0), r(N), p(N), Ap(N),is_bnd(N),g(N),b(N),ug(N)
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, A);
#else
	build_P1_CSRPattern(m, P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, A);
#endif
	//vol = M.sum();
	inited = false;
	iterate = 0;
	converged = false;
	lifted = false;
	
	compute_boundary_vertices(m,is_bnd);
	
	
	for (size_t i = 0; i < N; ++i) {

    double x = m.positions[i].x;
    double y = m.positions[i].y;

    // RHS
    b[i] =0;

    // non_homogenenus boundary ( affine)
    if (is_bnd[i]) {
		//g[i] = 0.0;
        g[i] = 1.0 + x + 2.0*y;
    } else {
        g[i] = 0.0;
    }
}
}

void PoissonSolver::clear_solution()
{
	lifted = false;
	#if USE_FEM_MATRIX
		build_P1_stiffness_matrix(m, A);
	#else
		build_P1_stiffness_matrix(m, P, A);
	#endif
	for (size_t i = 0; i < N; i++) {
		u[i] = 0;
	}
	init_cg();
	iterate = 0;
	converged = false;
}

//If i've BCs my problem is well posed, I don't need to set zero mean

/*void PoissonSolver::set_zero_mean(double *V)
{
	M.mvp(V, Ap.data);
	double s = blas_sum_in_place(Ap.data, N);
	for (size_t i = 0; i < N; ++i) {
		V[i] -= s / vol;
	}
}*/ 




void PoissonSolver::init_cg()
{
	double *F = f.data;//f evaluated in nodes
	double *U = u.data;//unknown
	double *R = r.data;//residue
	double *Pv = p.data;
	double *AP = Ap.data;//A*p
	lifted = false;
	// 0) we build A, we will modify it imposing boundary condition 
	#if !USE_FEM_MATRIX
	build_P1_stiffness_matrix(m, P, A);
	#endif
	
    
    double *B = b.data;// right-hand side values at vertices
    M.mvp(F, B);// B = M f

	//all I need to move on A
	uint32_t *rs = A.row_start;
	uint32_t *cc = A.col;
	double   *Ad = A.data.data;//It is A.data.data because A.data is a TArray<double> and I need a pointer to the data

	// 1) I set Ug[i] = g[i] for boundary vertices, (rememebering U[i]=0 for internal vertices)
	double *Ug = ug.data;   // new vector

	for (size_t i = 0; i < N; ++i) {

		if (is_bnd[i]) {
			Ug[i] = g[i];
		} else {
			Ug[i] = 0.0;
		}
    	U[i]  = 0.0;   // incognita correttiva iniziale
	}	

	// 2) B = B - A*Ug
    A.mvp(Ug, R);
    blas_axpy(-1.0, R, B, N);

	// 3) I set columns  of boundary vertices to 0, if not diagonal. 
	for (uint32_t i = 0; i < (uint32_t)N; ++i) {
		if (is_bnd[i]) continue;//skip internal vertices

		
		for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
			uint32_t j = cc[k];
			if (is_bnd[j]) Ad[k] = 0.0;
		}
		
	}
	// 3) I set rows  of boundary vertices equal to identity. 
	for (uint32_t i = 0; i < (uint32_t)N; ++i) {
		if (!is_bnd[i]) continue;//skip external vertices rows

		for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
			Ad[k] = 0.0;
		}
		
		for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
			
			if (cc[k] == i) {Ad[k] = 1.0;
			break;}
		}
		B[i] = 0;// sto imponenndo u0=0 sul bordo.
	}

// --- 5) initialize CG un A U = B

    // b2 = ||B||^2 (for relative error)
    b2 = blas_dot(B, B, N);

    
    blas_copy(B, R, N);   // perché U=0
  
    // P = R
    blas_copy(R, Pv, N);

    r2 = blas_dot(R, R, N);
    rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);

    inited = true;
}


void PoissonSolver::do_iterate(size_t max_iter, double tol)
{
	if (!inited) {
		init_cg();
	}

	double *U = u.data;
	double *R = r.data;
	double *Pv = p.data;
	double *AP = Ap.data;

	while (max_iter-- && rel_error > tol) {
		r2 = cg_iterate_once(A, U, R, Pv, AP, r2);
		iterate++;
		rel_error = sqrt(r2 / b2);
		if (rel_error <= tol) {
			converged = true;

			if (!lifted) {
				double *Ug = ug.data;
				for (size_t i = 0; i < N; ++i)
					U[i] += Ug[i];
				lifted = true;
			}
		}
	}

	if (rel_error <= tol)
		converged = true;
}
