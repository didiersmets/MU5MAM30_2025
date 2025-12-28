#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

static void load_cube_vertices(Vec3 *pos, size_t subdiv);
static void load_cube_indices(uint32_t *idx, size_t subdiv);

int load_cube(Mesh &m, size_t subdiv)
{
	/* Check subdiv is reasonable and return error if not */
	if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */) {
		return (-1);
	}

	size_t n = subdiv + 1;/*n is the number of vertices per face*/

	/* Reserve memory for vertices and indices */
	m.positions.resize(6 * POW2(n)); /*6 faces, each with n^2 vertices--> it duplicates por triplicates vertices which are respectively on edge or angle of cube*/
	m.indices.resize(36 * POW2(subdiv));/*6faces*3vertex per triangle*2 triagle per square*n^2 square per faces*/

	/* First build vertices as six unattached faces of n^2 vertices each */
	/* See below for implementation */
	load_cube_vertices(m.positions.data, subdiv);

	/* Build corresponding triangulation indices */
	/* See below for implementation */
	load_cube_indices(m.indices.data, subdiv);

	/* Finally attach faces between themselves */
	/* Implementation in src/duplicate_verts.cpp */
	remove_duplicate_vertices(m);

	return (0);
}

static void load_cube_vertices(Vec3 *pos, size_t subdiv)
{
	k=0;
	for (size_t f=0;f<6;f++){
		 for (size_t i=0;i<subdiv+1;i++){
			for (size_t j=0;j<subdiv+1;j++){
		
				if (f==0){/*face -X*/
					pos[(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j]={-1, -1 + 2* j/(double)subdiv, -1 + i*2/(double)subdiv};
				}
				else if (f==1){/*face -Y*/
					pos[(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j]={ -1 + 2* j/(double)subdiv,-1, -1 + i*2/(double)subdiv};
				}
				else if (f==2){/*face -Z*/
					pos[(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j]={ -1 + 2* j/(double)subdiv, -1 + i*2/(double)subdiv, -1};
				}
				else if (f==3){/*face +X*/
					pos[(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j]={ 1, -1 + j*2/(double)subdiv, -1 + 2* i/(double)subdiv};
				}
				else if (f==4){/*face +Y*/
					pos[(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j]={ -1 + 2* j/(double)subdiv, 1, -1 + i*2/(double)subdiv};
				}
				else if (f==5){/*face +Z*/
					pos[(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j]={ -1 + 2* j/(double)subdiv, -1 + i*2/(double)subdiv, 1};
				}

			}
		}	
	}
}

static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	/* Your implementation goes here */
	k=0;
	for (size_t f=0;f<6;f++){
		 for (size_t i=0;i<subdiv;i++){
			for (size_t j=0;j<subdiv;j++){
		
				/*ora mi trovo in un quadratino di una face */
				/*triangolo 1*/
				idx[k++]=(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j;
				idx[k++]=(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j+1;
				idx[k++]=(subdiv+1)*(subdiv+1)*f+(subdiv+1)*(i+1)+j+1;
				/*triangolo 2*/
				idx[k++]=(subdiv+1)*(subdiv+1)*f+(subdiv+1)*i+j;
				idx[k++]=(subdiv+1)*(subdiv+1)*f+(subdiv+1)*(1+i)+j;
				idx[k++]=(subdiv+1)*(subdiv+1)*f+(subdiv+1)*(1+i)+j+1;
				
			}
		}	
	}	
}
