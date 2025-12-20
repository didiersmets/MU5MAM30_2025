#include <stdio.h>
#include <fstream>

#include "mesh.h"

using namespace std;

void save_to_obj(Mesh &m, string file_name){
	ofstream ofs {file_name};
	for (size_t i = 0; i<m.vertex_count(); i++){
		ofs << "v " << m.positions[i].x << " " << m.positions[i].y << " " << m.positions[i].z << "\n";
	}
	// warning vertices are indexed starting from 1 (https://en.wikipedia.org/wiki/Wavefront_.obj_file#Face_elements)
	for (size_t i = 0; i< m.index_count(); i+=3){
		ofs << "f " << (m.indices[i]+1) << " " << (m.indices[i+1]+1) << " " << (m.indices[i+2]+1) << "\n";
	}
}
