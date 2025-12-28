#include <stdio.h>
#include "mesh.h"
#include "sphere.h"
#include "P1.h"

int main(){
	Mesh m;
	load_sphere(m,2);
	CSRPattern p;
	build_P1_CSRPattern(m,p);
	spy(p,p.rows,"CSRPattern.png");
	CSRMatrix M(p);
	CSRMatrix S(p);
	build_P1_mass_matrix(m,M);
	build_P1_stiffness_matrix(m,S);

}
