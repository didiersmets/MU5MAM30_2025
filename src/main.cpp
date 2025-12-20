#include <stdio.h>
#include "mesh.h"
#include "sphere.h"

int main(){
	Mesh m;
	load_sphere(m,4);
	save_to_obj(m,"sphere.obj");
}
