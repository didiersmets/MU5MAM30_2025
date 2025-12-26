#include <stdio.h>
#include "mesh.h"
#include "sphere.h"
#include "P1.h"

int main(){
	Mesh m;
	load_sphere(m,2);
	CSRPattern p;
	build_P1_CSRPattern(m,p);
	CSRMatrix mat(p);
	spy(p,p.rows,"CSRPattern.png");

}
