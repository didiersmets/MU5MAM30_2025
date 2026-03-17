#include <stdio.h>
#include "mesh.h"
#include "sphere.h"
#include "P1.h"
#include "conjugate_gradient.h"
#include "sparse_tri_solve.h"

int main(){
	Mesh m;
	load_sphere(m,0);
	CSRPattern pattern {0};
	build_P1_CSRPattern(m,pattern);
	size_t n = pattern.rows; // square matrix
	spy(pattern,pattern.rows,"CSRPattern.png");
	CSRMatrix A(pattern,0);

	// we want to solve Mx = b
	dump(A,"matrix0.txt");
	build_P1_mass_matrix(m,A);
	dump(A,"matrix1.txt");
	build_P1_stiffness_matrix(m,A);
	dump(A,"matrix2.txt");



	TArray<double> b (n,0.0);
	b[0] = 1.0;
	TArray<double> x (n,0.0);
	TArray<double> r (n);
	TArray<double> p (n);
	TArray<double> Ap (n);
	double rel_error = 0;
	size_t iter = conjugate_gradient_solve(A,b.data,x.data,r.data,p.data,Ap.data,&rel_error,1e-6,n,true);
	printf("Solved in %ld iteration (relative error %lf)\n",iter,rel_error);
}
