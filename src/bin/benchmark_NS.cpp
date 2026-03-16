#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <chrono>

#include "cube.h"
#include "sphere.h"
#include "mesh.h"
#include "mesh_io.h"
#include "navier_stokes.h"

static void syntax(const char* prg)
{
    printf("Syntax:\n");
    printf("  %s <cube|cube2|sphere|obj> <n_or_filename> [steps] [dt] [nu]\n", prg);
}

static int load_mesh_from_args(Mesh& mesh, const char* kind, const char* arg)
{
    if (strcmp(kind, "cube") == 0)  return load_cube(mesh, atoi(arg));
    if (strcmp(kind, "cube2") == 0) return load_cube_nested_dissect(mesh, atoi(arg));
    if (strcmp(kind, "sphere") == 0) return load_sphere(mesh, atoi(arg));
    if (strcmp(kind, "obj") == 0)   return load_obj(arg, mesh);
    return -1;
}

static void init_vorticity(NavierStokesSolver& solver)
{
    for (size_t i = 0; i < solver.N; ++i) {
        double x = solver.m.positions[i].x;
        double y = solver.m.positions[i].y;
        double z = solver.m.positions[i].z;

        solver.omega[i] =
            0.8 * exp(-6.0 * (x * x + y * y))
            + 0.3 * sin(3.0 * x)
            - 0.2 * cos(4.0 * z);
    }

    solver.set_zero_mean(solver.omega.data);
    memset(solver.psi.data, 0, solver.N * sizeof(double));
    solver.t = 0.0;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        syntax(argv[0]);
        return EXIT_FAILURE;
    }

    const char* mesh_kind = argv[1];
    const char* mesh_arg  = argv[2];

    int steps   = (argc > 3) ? atoi(argv[3]) : 200;
    double dt   = (argc > 4) ? atof(argv[4]) : 0.002;
    double nu   = (argc > 5) ? atof(argv[5]) : 2e-4;

    Mesh mesh;
    if (load_mesh_from_args(mesh, mesh_kind, mesh_arg) != 0) {
        printf("Error: could not load mesh.\n");
        return EXIT_FAILURE;
    }

    NavierStokesSolver solver(mesh);
    init_vorticity(solver);

    double total = 0.0;
    double first_step_time = 0.0;
    double rest_steps_total = 0.0;
    double min_step = 1e100;
    double max_step = 0.0;

    for (int k = 0; k < steps; ++k) {
        auto tk0 = std::chrono::high_resolution_clock::now();

        solver.time_step(dt, nu);

        auto tk1 = std::chrono::high_resolution_clock::now();
        double step_time = std::chrono::duration<double>(tk1 - tk0).count();

        total += step_time;

        if (k == 0) {
            first_step_time = step_time;
        } else {
            rest_steps_total += step_time;
        }

        if (step_time < min_step) min_step = step_time;
        if (step_time > max_step) max_step = step_time;
    }

    double avg_all  = total / (double)steps;
    double avg_rest = (steps > 1) ? rest_steps_total / (double)(steps - 1) : 0.0;

    double minv = solver.omega[0];
    double maxv = solver.omega[0];
    double norm2 = 0.0;
    for (size_t i = 0; i < solver.N; ++i) {
        double v = solver.omega[i];
        if (v < minv) minv = v;
        if (v > maxv) maxv = v;
        norm2 += v * v;
    }

    printf("=====================================\n");
#if NS_USE_DIRECT_SOLVER
    printf("mode         = Cholesky\n");
#else
    printf("mode         = CG\n");
#endif
    printf("mesh         = %s %s\n", mesh_kind, mesh_arg);
    printf("DoF          = %zu\n", solver.N);
    printf("steps        = %d\n", steps);
    printf("dt           = %.6e\n", dt);
    printf("nu           = %.6e\n", nu);
    printf("total time   = %.12e s\n", total);
    printf("first step   = %.12e s\n", first_step_time);
    printf("avg/step     = %.12e s\n", avg_all);
    printf("avg/rest     = %.12e s\n", avg_rest);
    printf("min step     = %.12e s\n", min_step);
    printf("max step     = %.12e s\n", max_step);
    printf("final t      = %.12e\n", solver.t);
    printf("omega min    = %.12e\n", minv);
    printf("omega max    = %.12e\n", maxv);
    printf("||omega||_2  = %.12e\n", sqrt(norm2));
    printf("=====================================\n");

    return EXIT_SUCCESS;
}