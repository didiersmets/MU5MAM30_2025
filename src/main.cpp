#include <stdio.h>
#include "mesh/mesh.h"
#include "mesh/sphere.h"

int main(){
	Mesh m;
	load_sphere(m,0);
	save_to_obj(m,"test_obj.obj");
	printf("Hello world\n");
}
