#include <stdio.h>
#include "mesh.h"
#include "sphere.h"
#include "P1.h"

int main(){
	Mesh m;
	load_sphere(m,4);
	CSRPattern p;
	build_P1_CSRPattern(m,p);
	CSRMatrix mat(p);

}
