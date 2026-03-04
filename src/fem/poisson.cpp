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
    : m(m), N(m.vertex_count()), f(N), u(N, 0.0), r(N), p(N), Ap(N)
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
	compute_boundary_vertices(m,is_bnd);
	g.resize(N);
	b.resize(N);
	for (size_t i = 0; i < N; ++i) {

    double x = veritces[i].x;
    double y = vertices[i].y;

    // RHS
    b[i] =0;

    // bordo non omogeneo esempio affine
    if (is_boundary[i]) {
		g[i] = 0.0;
       // g[i] = 1.0 + x + 2.0*y;
    } else {
        g[i] = 0.0;
    }
}
}

void PoissonSolver::clear_solution()
{
	build_P1_stiffness_matrix(m, P, A);
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
	double *F = f.data;
	double *U = u.data;
	double *R = r.data;
	double *Pv = p.data;
	double *AP = Ap.data;
	
	/* Fix up F and U for zero mean */
	/*set_zero_mean(F);
	set_zero_mean(U);*/

	//1--I set buondary conditions
	for (size_t i = 0; i < N; ++i) {
        if (is_bnd[i]) {
            U[i] = g[i];
        }
    }
	//2--I compute and b2
	 if (b.size != N) {
        b.resize(N);
    }
    double *B = b.data;
    M.mvp(F, B);//B = Mf

// Puntatori CSR per A
    uint32_t *rs = A.row_start;//tells me where each row starts in the col and data arrays
    uint32_t *cc = A.col;//tells me the column index of each nonzero entry in A
    double   *Ad = A.data.data;//nonnzero values of A

    // --- 3) Correzione RHS: B_i -= sum_{j boundary} A_ij * g_j
    // (si applica bene per tutti; per i di bordo poi sovrascriviamo comunque)
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        double corr = 0.0;
		//nel secondo for sta facendo il profotto tra A e ug (dove ug è zero se non sono di bordo e g[i] se sono di bordo) e lo sta sottraendo a B[i]
        for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {//evaluate elements in row i
            uint32_t j = cc[k];//column index of the current element
            if (is_bnd[j]) {
                corr += Ad[k] * g[j];
            }
        }
        B[i] -= corr;
    }
	// --- 3b) Per mantenere A simmetrica per CG:
	// per ogni riga interna i, azzera A(i,j) se j è di bordo.
	for (uint32_t i = 0; i < (uint32_t)N; ++i) {
   	 if (is_bnd[i]) continue; // riga di bordo la sistemo dopo

    	for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
        uint32_t j = cc[k];
        if (is_bnd[j]) {
            Ad[k] = 0.0;
        }
    }
}
    // --- 4) Impone Dirichlet nella matrice e nel RHS:
    // per ogni i di bordo: riga i = e_i^T, B_i = g_i
    for (uint32_t i = 0; i < (uint32_t)N; ++i) {
        if (!is_bnd[i]) continue;

        // azzera riga i
        for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
            Ad[k] = 0.0;
        }

        // metti diagonale A(i,i)=1 (assumiamo che (i,i) sia nel pattern, cosa vera per stiffness P1)
        bool found_diag = false;
        for (uint32_t k = rs[i]; k < rs[i + 1]; ++k) {
            if (cc[k] == i) {
                Ad[k] = 1.0;
                found_diag = true;
                break;
            }
        }

        B[i] = g[i];// sta sovrascrivendo solo il valore di B[i] per i di bordo, mettendolo uguale a g[i]
    }


	// --- 5) Ora inizializza CG su A U = B

    // b2 = ||B||^2 (per errore relativo)
    b2 = blas_dot(B, B, N);

    // R = B - A U
    A.mvp(U, R);//R = AU;
    blas_axpy(-1.0, R, B, N);   // B = B - AU (uso B come temporaneo)
    blas_copy(B, R, N);         // R = B - AU

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
	}

	if (rel_error <= tol)
		converged = true;
}
