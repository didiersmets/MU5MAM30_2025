#pragma once

#include <stddef.h>
#include <vector>
#include "matrix.h"

// 单步投影共轭梯度迭代
double projected_cg_iterate_once(const Matrix &A, double *__restrict x,
                                 double *__restrict r, double *__restrict p,
                                 double *__restrict Ap, double r2,
                                 bool has_boundary, 
                                 const std::vector<bool> *is_boundary);

// 投影共轭梯度求解器主控函数
size_t projected_conjugate_gradient_solve(const Matrix &A, const double *__restrict b,
                                          double *__restrict x, double *__restrict r,
                                          double *__restrict p, double *__restrict Ap,
                                          double *rel_error, double tol, int max_iter,
                                          bool inited,
                                          bool has_boundary, 
                                          const std::vector<bool> *is_boundary);