#pragma once

#include <functional>

#include "array.h"
#include "my_sparse_matrix.h"
#include "my_mesh.h"
#include "vec2.h"

// Vec2 equal to TVec2<double>
using Vec2d = TVec2<double>;

struct MyMinimalGraphSolver {
    MyMinimalGraphSolver(const MyMesh &m, std::function<double(const Vec2d&)> func);
    const MyMesh &m;
    size_t N;   
    size_t N_b; 

    TArray<double> u; 
    TArray<double> uold;
    TArray<double> du; 
    TArray<double> q; 

    TArray<double> f; 
    TArray<double> b; 
    
    TArray<double> r;  
    TArray<double> p; 
    TArray<double> Ap; 

    CSRPattern P;
    CSRMatrix S, S_modified;

    bool inited; 
    size_t iterate_N, iterate_P; 
    TArray<double> residual_N, residual_P; 
    bool converged;

    static const size_t iter_max = 500;
    static constexpr double toll = 1e-12;
    const double tolCG = 1e-16;
    
    void clear_solution(bool Newton);
    void do_iterate_Newton(size_t max_iter = iter_max, double tol = toll, const double min_alpha = 1.0, const double c = 1e-4, const double rho = 0.5);
    void do_iterate_Picardi(size_t max_iter = iter_max, double tol = toll);
    double compute_denominator(TArray<double> &den, const TArray<double> &u); //returns the area
};