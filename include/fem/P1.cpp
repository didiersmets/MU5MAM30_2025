

#include <stdio.h>
#include <string.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "mass.h"
#include "stiffness.h"
#include "sys_utils.h"
#include "mesh.h"
#include "sparse_matrix.h"


void build_P1_mass_matrix(const Mesh &m, FEMatrix &M)
{
   
}

void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S)
{
}

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
}

void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
}

int main(){
   
    
Mesh mesh;
load_cube(mesh, 0);  
remove_duplicate_vertices(mesh);

FEMatrix M;
build_P1_mass_matrix(mesh, M);

// Check properties:
// 1. Diagonal entries should be positive
// 2. Off-diagonal entries should be positive (for mass matrix)
// 3. Matrix sum should equal total mesh volume / 3 (roughly)
std::cout << "Matrix sum: " << M.sum() << std::endl;


}