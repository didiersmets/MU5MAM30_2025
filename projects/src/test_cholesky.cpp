#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip> // for std::setprecision

#include "mesh.h"
#include "mesh_bounds.h"
#include "mesh_gpu.h"
#include "mesh_io.h"
#include "cube.h"
#include "symbolic.h"
#include "sparse_matrix.h"
#include "array.h"
#include "P1.h"
#include "cholesky.h"
#include "poisson.h"
#include "conjugate_gradient.h"

#include "tiny_expr/tinyexpr.h"

/* RHS expression of the PDE */
char rhs_expression[128] =
	"cos(35 * y * sin(27 + 13 * x^2 + 19 * z^2 - 13 * x * z))";
bool rhs_show_error = false;
double rhs_x, rhs_y, rhs_z, rhs_p, rhs_t, rhs_r;
te_variable rhs_vars[] = { { "x", &rhs_x },	{ "y", &rhs_y },
			   { "z", &rhs_z },	{ "phi", &rhs_p },
			   { "theta", &rhs_t }, { "rand", &rhs_r } };
te_expr *te_rhs = NULL;

bool new_rhs_chol(double* f, Mesh &m)
{
	te_expr *test = te_compile(rhs_expression, rhs_vars,
				   sizeof(rhs_vars) / sizeof(rhs_vars[0]),
				   NULL);
	if (!test)
		return false;

	te_free(te_rhs);
	te_rhs = test;
	for (size_t i = 0; i < m.vertex_count(); ++i) {
		rhs_x = m.positions[i].x;
		rhs_y = m.positions[i].y;
		rhs_z = m.positions[i].z;
		rhs_p = atan2(rhs_y, rhs_x);
		rhs_t = atan2(sqrt(rhs_x * rhs_x + rhs_y * rhs_y), rhs_z);
		rhs_r = (double)rand() / RAND_MAX;
		f[i] = te_eval(te_rhs);
	}

	return true;
}

//convert matrix in CSRFormat to dense
void CSR_to_dense(CSRMatrix &M, TArray<double> &M_dense){
  size_t n = M.rows;
  assert(M_dense.size == n*n);

  //first fill with zeros
  for (size_t i = 0; i < n*n; ++i) {
    M_dense[i] = 0;}

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = M.row_start[i]; j < M.row_start[i+1]; j++) {
        M_dense[i*n + M.col[j]] = M.data[j];
        if(M.symmetric){
            M_dense[M.col[j]*n + i] = M.data[j];
        }   
    }
    }
    
}

// Print a dense matrix stored in a vector of size n*n
void print_dense_matrix(const TArray<double> &M, size_t n, int precision = 4) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << std::setw(10) << std::fixed 
                      << std::setprecision(precision) 
                      << M[i * n + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]){
    bool do_checks = false;
    size_t subdiv = 1;

    // Parse command line arguments
    //Command line argument parsing from ChatGPT
    for(int i = 1; i < argc; i++){
        std::string arg = argv[i];
        if(arg == "--check" || arg == "-c") {
            do_checks = true;
        } else if((arg == "--subdiv" || arg == "-s") && i + 1 < argc) {
            subdiv = static_cast<size_t>(std::atoi(argv[i+1]));
            i++; // skip next arg
        }
    }
    //Load mesh
    Mesh mesh;
    load_cube(mesh, subdiv);

    CSRPattern P, L_pattern;
    build_P1_CSRPattern(mesh, P);


    TArray<uint32_t> parent(P.rows);
    construct_etree(P, parent);
    construct_L_sparsity_pattern(P, L_pattern, parent);

    CSRMatrix M, L, A, M_check;

    build_P1_mass_matrix(mesh, P, M);


    cholesky_factorization(M, L_pattern, L);
    size_t n = P.rows;
    TArray<double> x(P.rows, 0.0),  f(P.rows, 0.0), f_res(P.rows, 0.0);

     if(do_checks){
    // Compute LL^T
    
    TArray<double> M_dense(n*n), L_dense(n*n), LLt(n*n);
    CSR_to_dense(M, M_dense);
    CSR_to_dense(L, L_dense);
    for(size_t i   = 0; i < n; i++){
        for(size_t j = 0; j < n; j ++){
            LLt[i*n + j] = 0;
            for(size_t k = 0; k <n; k ++){
                LLt[i*n + j] += L_dense[i*n +k] * L_dense[j*n +k];
            }
        }
    }

    //compute the frobenious norm
    double err = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double d = M_dense[i*n +j ] - LLt[i*n +j ];
            err += d * d;
        }
    }
    err = sqrt(err);
    std::cout << "Error " << err << std::endl;
   
    new_rhs_chol(f.data, mesh);
    solve_cholesky(L, f.data, x.data);

    M.mvp(x.data, f_res.data); 
    err = 0;
    for (size_t j = 0; j < n; ++j) {
            double d = f[j ] - f_res[j ];
            err += d * d;
    }

        std::cout << "Error " << err << std::endl;
    }

    // Time solving if correctness checks are disabled
    if(!do_checks){
        std::cout << "subdiv: " << subdiv << std::endl;
        std::cout << "n: " << n << std::endl;

        auto t0 = std::chrono::high_resolution_clock::now();
        
        construct_etree(P, parent);
        construct_L_sparsity_pattern(P, L_pattern, parent);
        cholesky_factorization(M, L_pattern, L);
        auto t1 = std::chrono::high_resolution_clock::now();
        solve_cholesky(L, f.data, x.data);
       
        std::chrono::duration<double> elapsed = t1 - t0;
        std::cout << "chol_time: " << elapsed.count() << " seconds\n";
        Mesh mesh_disect;
        load_cube_nested_dissect(mesh_disect, subdiv);
    
        //same experiments but using nested disection
        CSRPattern P_D, LD_pattern;
        build_P1_CSRPattern(mesh_disect, P_D);
        TArray<uint32_t> parent_D(P_D.rows);
        build_P1_mass_matrix(mesh_disect, P_D, M);


        t0 = std::chrono::high_resolution_clock::now();
        construct_etree(P_D, parent_D);
        construct_L_sparsity_pattern(P_D, LD_pattern, parent_D);
        cholesky_factorization(M, LD_pattern, L);
        std::cout << "chol_nested_time: " << elapsed.count() << " seconds\n";
        solve_cholesky(L, f.data, x.data);
         t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
       
    } 

    return 0;



}