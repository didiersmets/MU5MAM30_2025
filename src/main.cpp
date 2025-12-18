#include <stdio.h>
#include "mesh.h"
#include "sphere.h"

int main(){
	Mesh m;
	load_sphere(m,0);
	save_to_obj(m,"test_obj.obj");
	printf("Hello world\n");
}
