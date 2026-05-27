#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <functional>

#include "my_minimal_graph.h"

#include "my_P1.h"
#include "tiny_blas.h"
#include "my_conjugate_gradient.h"

//#include "export_mesh.h"

#define HUGE 1e30
using Vec2d = TVec2<double>;

MyMinimalGraphSolver::MyMinimalGraphSolver(const MyMesh &m, std::function<double(const Vec2d&)> func)
    : m(m), N(m.vtx_count), N_b(m.boundary_count), u(N), uold(N), du(N), q(m.tri_count), f(N),
        b(N, 0.0), r(N), p(N), Ap(N), inited(false), iterate_N(0), iterate_P(0), residual_N(iter_max,0.0), 
        residual_P(iter_max,0.0), converged(false){
    
    build_P1_CSRPattern(m, P); 
    
    // assign values to f using the provided function func, which takes a Vec2d as input and returns a double
    for (size_t i = 0; i < N; i++) {
        f[i] = func(m.vertices[i]); // assign value of vertex i with 2D coordinates m.vertices[i] to f[i]
    }

}

void MyMinimalGraphSolver::clear_solution(bool Newton)
{
    memset(b.data, 0, N * sizeof(double)); // clear right hand side
    memset(du.data, 0, N * sizeof(double)); // clear solution increment

    // set boundary values b to f
    for (size_t i = 0; i < N_b; i++) {
        b[m.boundary[i]] = f[m.boundary[i]];
    }
    memcpy(u.data, b.data, N * sizeof(double)); // set initial solution u to boundary values, this will be used as the initial guess for the solver
    memcpy(uold.data, b.data, N * sizeof(double)); // set uold to boundary values, this will be used to compute the residual for the stopping criterion
    converged = false;
    inited = true;
    if(Newton)
        iterate_N = 0;
    else
        iterate_P = 0;
}

double MyMinimalGraphSolver::compute_denominator(TArray<double> &den, const TArray<double> &u){

    /*
    compute the denominator of the energy functional
    */
   
    uint32_t a, b, c;
    Vec2d A, B, C, AB, AC;
    double S_loc[6], area = 0.0, tri_area, ABAB, ACAC, ABAC, mult;
    for(size_t t = 0; t < m.tri_count; t++){
        a = m.triangles[t].x;
        b = m.triangles[t].y;
        c = m.triangles[t].z;
        A = m.vertices[a];
		B = m.vertices[b];
		C = m.vertices[c];
		AB = {(double)B[0] - (double)A[0],
			    (double)B[1] - (double)A[1]};
		AC = {(double)C[0] - (double)A[0],
			    (double)C[1] - (double)A[1]};

        ABAB = dot(AB,AB);
	    ACAC = dot(AC,AC);
	    ABAC = dot(AB, AC);
        tri_area = 0.5 * sqrt(ABAB * ACAC - ABAC * ABAC);
	    mult = 0.25 / (tri_area);
	    ABAB *= mult;
	    ACAC *= mult;
	    ABAC *= mult;

	    S_loc[0] = ACAC - 2 * ABAC + ABAB;
	    S_loc[1] = ACAC;
	    S_loc[2] = ABAB;
        S_loc[3] = ABAC - ACAC;
	    S_loc[4] = -ABAC;
	    S_loc[5] = ABAC - ABAB;

        den[t] = 1.0 / sqrt(1 + u[a]*u[a]*S_loc[0] + u[b]*u[b]*S_loc[1] + u[c]*u[c]*S_loc[2] +
                            2*(u[a]*u[b]*S_loc[3] + u[b]*u[c]*S_loc[4] + u[c]*u[a]*S_loc[5]));
        area += tri_area/den[t];
    }
    return area;
}



void MyMinimalGraphSolver::do_iterate_Newton(size_t max_iter, double tol, const double min_alpha, const double c, const double rho)
{

    // we have already cleared before calling this function, so we can start iterating directly
    if (!inited) {
        clear_solution(true);
    }
    TArray<double> u_tmp(N, 0.0);
    int iterCG;
    double error2 = 0.0, errorCG = 0.0, area;
    double alpha = 1.0, energy_tmp = 0.0;
    bool flag = true;
    double prev_error = 1e30;  // Track previous error for stagnation detection
    int stagnant_count = 0;    // Count consecutive iterations without improvement

    // set initial guess to boundary values, this will be used as the initial guess for the solver
    // we have already set u to boundary values in clear_solution, so we can skip this step
    
    //memcpy(u.data,b.data ,N * sizeof(double)); 
    
    
    //export_mesh_to_csv(m, u, "start.csv");

    area = compute_denominator(q, u); // computes denominator q for the initial solution u
    if (iterate_N == 0) {
        printf("Starting Newton solver.... \n");
        printf("%-10s %-15s %-15s %-15s %-15s\n", "Iter", "ErrorNewton", "IterCG", "ErrorCG", "Area");
        printf("%-10s %-15s %-15s %-15s %-15g \n", "-", "-", "-", "-", area);
    }
    size_t target_iter = iterate_N + max_iter;
    while (iterate_N < target_iter){
        
        build_P1_stiffness_matrix_NS(m, P, S_modified, q.data, u.data); // build the jacobian stiffness matrix S_modified for the current solution u, using the denominator q computed from u
        build_P1_rhs_NS(m, q.data, u.data, b); // build the right hand side b for the current solution u, using the denominator q computed from u

        for (size_t i = 0; i < N_b; i++){
            S_modified(m.boundary[i],m.boundary[i]) = HUGE; // set boundary values to a large number in the stiffness matrix, this will enforce the solution to be equal to the boundary values at the boundary vertices
            b[m.boundary[i]] = 0; // set rhs boundary values to 0 in the right hand side, such that solution increment at boundary will be 0
        }

        // solve the linear system S_modified * du = b 
        // du is the solution increment, we will update the solution u by adding alpha * du, where alpha is a step size determined by a line search
        iterCG = conjugate_gradient_solve(S_modified, b.data, du.data, r.data, p.data, Ap.data , &errorCG, tolCG, 10000, false); 

        area = compute_denominator(q,u);
        error2 = blas_dot(du.data, du.data, N); // compute the squared norm of the solution increment du, this will be used as the error for the stopping criterion

        flag = true;
        alpha = 1.0;

        while(flag){
            // update each entry of u_tmp with the current solution u plus alpha times the solution increment du
            for (size_t i = 0; i < N; i++){
                // u_tmp[i] = u[i] + alpha*du[i]; // update with alpha times the solution increment, alpha could be removed, this option here seems to converge faster
                u_tmp[i] = u[i] + du[i]; // update with alpha times the solution increment, alpha could be removed

            }
            energy_tmp = compute_denominator(q,u_tmp);
            // this is the condition for alpha of the Armijo backtracking line search
            // maybe we can remove that
            if (energy_tmp < area - c*alpha*error2 || alpha <= min_alpha){
                flag = false;
            }
            else{
                alpha *= rho;
            }
        }
        
        for (size_t i = 0; i < N; i++){
            u[i] += alpha*du[i];
        }

        
        error2 = sqrt(error2);
        
        printf("%-10ld %-15g %-15d %-15g %-15g\n", iterate_N, error2, iterCG, errorCG, area);

        // Store residual only if within array bounds
        if (iterate_N < residual_N.size) {
            residual_N[iterate_N] = error2;
        }
        iterate_N++;
        
        // Check for convergence
        if (error2 < tol) {
            converged = true;
            printf("Converged after %ld iterations.\n", iterate_N);
            break;
        }
        
        // Stagnation detection: if error hasn't improved significantly for 10 iterations, stop
        if (error2 > prev_error * 0.9999) {  // Less than 0.01% improvement
            stagnant_count++;
            if (stagnant_count > 10) {
                printf("Stagnation detected after %ld iterations (error not decreasing).\n", iterate_N);
                converged = false;
                break;
            }
        } else {
            stagnant_count = 0;  // Reset if we make progress
        }
        
        prev_error = error2;
        
        // Prevent running beyond the maximum iteration limit (avoid accessing residual array out of bounds)
        if (iterate_N >= iter_max) {
            printf("Reached maximum iteration limit (%zu).\n", iter_max);
            break;
        }
        
    }
    if (!converged && error2 <= tol) {
        converged = true;
        printf("Converged after %ld iterations.\n", iterate_N);
    }
    else if (!converged) {
        printf("Did not converge after %ld iterations.\n", iterate_N);
    }

    //export_mesh_to_csv(m, u, "solutionNewton.csv");
    //save_residual_to_csv(residual_N, "residualNewton.csv");
}

/*
void::MyMinimalGraphSolver::do_iterate_Picardi(size_t max_iter, double tol)
{   
    clear_solution(false);
    int iterCG;
    double error2 = 0.0, errorCG = 0.0, area;

    export_mesh_to_csv(m, u, "start.csv");
   
    for (size_t i = 0; i < N_b; i++){
        b[m.boundary[i]] *= HUGE;
    }

    area = compute_denominator(q,u);

    printf("Starting Picardi Solver.... \n");
    printf("%-10s %-15s %-15s %-15s %-15s\n", "Iter", "ErrorPicardi", "IterCG", "ErrorCG", "Area");
    printf("%-10s %-15s %-15s %-15s %-15g \n", "-", "-", "-", "-", area);
    while (iterate_P < max_iter){
        
        memcpy(uold.data, u.data, N * sizeof(double));

        build_P1_stiffness_matrix(m, P, S, true, q.data);

        for (size_t i = 0; i < N_b; i++){
            S(m.boundary[i],m.boundary[i]) = HUGE;           
        }

        iterCG = conjugate_gradient_solve(S, b.data, u.data, r.data, p.data, Ap.data , &errorCG, tolCG, 10000, false);

        error2 = 0.0;

        for(size_t i = 0; i < N; i++)
            error2 += (u[i] - uold[i])*(u[i] - uold[i]);
        
        error2 = sqrt(error2);

        area = compute_denominator(q,u);

        printf("%-10ld %-15g %-15d %-15g %-15g\n", iterate_P, error2, iterCG, errorCG, area);

        residual_P[iterate_P] = error2;

        iterate_P++;
        if (error2 < tol)
            break;
        
    }
    residual_P.resize(iterate_P);
    if (error2 <= tol) {
        converged = true;
        printf("Converged after %ld iterations.\n", iterate_P);
    }

    else {
        printf("Did not converge after %ld iterations.\n", iterate_P);
    }
    export_mesh_to_csv(m, u, "solutionPicardi.csv");
    save_residual_to_csv(residual_P, "residualPicardi.csv");
}
*/