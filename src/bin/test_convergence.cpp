#include "mesh.h"
#include "fem/poisson.h"
#include "fem/fem_type.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <utility>
#include <algorithm>


// ------------------------------
// MALLA [0,1]^2
// ------------------------------
void build_unit_square_mesh(Mesh &mesh, int N)
{
    mesh.positions.clear();
    mesh.indices.clear();
    mesh.attr.clear();

    mesh.positions.resize((N + 1) * (N + 1));

    double h = 1.0 / N;

    for (int j = 0; j <= N; ++j) {
        for (int i = 0; i <= N; ++i) {
            size_t id = j * (N + 1) + i;
            mesh.positions[id] = Vec3(i * h, j * h, 0.0);
        }
    }

    mesh.indices.reserve(6 * N * N);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {

            int v0 = j * (N + 1) + i;
            int v1 = v0 + 1;
            int v2 = v0 + (N + 1);
            int v3 = v2 + 1;

            mesh.indices.push_back(v0);
            mesh.indices.push_back(v1);
            mesh.indices.push_back(v2);

            mesh.indices.push_back(v1);
            mesh.indices.push_back(v3);
            mesh.indices.push_back(v2);
        }
    }

    mesh.attr.resize(mesh.vertex_count());
    for (size_t i = 0; i < mesh.vertex_count(); ++i)
        mesh.attr[i] = 0.0f;
}

// ------------------------------
// SOLUCIÓN EXACTA
// ------------------------------

double exact_u(double x, double y)
{
    return sin(2*M_PI*x) * sin(2*M_PI*y);
}

double exact_f(double x, double y)
{
    return 8*M_PI*M_PI * sin(2*M_PI*x) * sin(2*M_PI*y);
}




// ------------------------------
// calculos de gradientes para P1 y P2
// ------------------------------

void exact_grad_u(double x, double y, double &ux, double &uy)
{
    ux = 2*M_PI * cos(2*M_PI*x) * sin(2*M_PI*y);
    uy = 2*M_PI * sin(2*M_PI*x) * cos(2*M_PI*y);
}


void grad_u_P1(const Mesh &mesh, const TArray<double> &u,
               int a, int b, int c,
               double &ux, double &uy)
{
    Vec3 A = mesh.positions[a];
    Vec3 B = mesh.positions[b];
    Vec3 C = mesh.positions[c];

    double x1=A.x, y1=A.y;
    double x2=B.x, y2=B.y;
    double x3=C.x, y3=C.y;

    double det = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);

    double beta1  = (y2 - y3) / det;
    double gamma1 = (x3 - x2) / det;

    double beta2  = (y3 - y1) / det;
    double gamma2 = (x1 - x3) / det;

    double beta3  = (y1 - y2) / det;
    double gamma3 = (x2 - x1) / det;

    ux = beta1*u[a] + beta2*u[b] + beta3*u[c];
    uy = gamma1*u[a] + gamma2*u[b] + gamma3*u[c];
}

// Devuelve lambda1, lambda2, lambda3 y gradientes de las baricéntricas
void barycentric_and_gradients(
    const Vec3 &A, const Vec3 &B, const Vec3 &C,
    double x, double y,
    double &l1, double &l2, double &l3,
    double grad_l1[2], double grad_l2[2], double grad_l3[2])
{
    double x1=A.x, y1=A.y;
    double x2=B.x, y2=B.y;
    double x3=C.x, y3=C.y;

    double detT = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);

    // Baricéntricas
    l1 = ((x2-x)*(y3-y) - (x3-x)*(y2-y)) / detT;
    l2 = ((x3-x)*(y1-y) - (x1-x)*(y3-y)) / detT;
    l3 = 1.0 - l1 - l2;

    // Gradientes de las baricéntricas
    grad_l1[0] = (y2 - y3) / detT;
    grad_l1[1] = (x3 - x2) / detT;

    grad_l2[0] = (y3 - y1) / detT;
    grad_l2[1] = (x1 - x3) / detT;

    grad_l3[0] = (y1 - y2) / detT;
    grad_l3[1] = (x2 - x1) / detT;
}

void eval_P2(
    const Mesh &mesh,
    const TArray<double> &u,
    int a, int b, int c, int ab, int bc, int ca,
    double x, double y,
    double &u_val,
    double &ux, double &uy)
{
    Vec3 A = mesh.positions[a];
    Vec3 B = mesh.positions[b];
    Vec3 C = mesh.positions[c];

    double l1, l2, l3;
    double g1[2], g2[2], g3[2];

    barycentric_and_gradients(A,B,C, x,y, l1,l2,l3, g1,g2,g3);

    // Funciones de forma P2
    double phi[6];
    phi[0] = l1*(2*l1 - 1);
    phi[1] = l2*(2*l2 - 1);
    phi[2] = l3*(2*l3 - 1);
    phi[3] = 4*l1*l2;
    phi[4] = 4*l2*l3;
    phi[5] = 4*l3*l1;

    // Gradientes de las funciones de forma
    double dphi[6][2];

    dphi[0][0] = (4*l1 - 1)*g1[0];
    dphi[0][1] = (4*l1 - 1)*g1[1];

    dphi[1][0] = (4*l2 - 1)*g2[0];
    dphi[1][1] = (4*l2 - 1)*g2[1];

    dphi[2][0] = (4*l3 - 1)*g3[0];
    dphi[2][1] = (4*l3 - 1)*g3[1];

    dphi[3][0] = 4*(l1*g2[0] + l2*g1[0]);
    dphi[3][1] = 4*(l1*g2[1] + l2*g1[1]);

    dphi[4][0] = 4*(l2*g3[0] + l3*g2[0]);
    dphi[4][1] = 4*(l2*g3[1] + l3*g2[1]);

    dphi[5][0] = 4*(l3*g1[0] + l1*g3[0]);
    dphi[5][1] = 4*(l3*g1[1] + l1*g3[1]);

    // Evaluar u_h
    u_val =
        phi[0]*u[a] +
        phi[1]*u[b] +
        phi[2]*u[c] +
        phi[3]*u[ab] +
        phi[4]*u[bc] +
        phi[5]*u[ca];

    // Evaluar gradiente
    ux =
        dphi[0][0]*u[a] +
        dphi[1][0]*u[b] +
        dphi[2][0]*u[c] +
        dphi[3][0]*u[ab] +
        dphi[4][0]*u[bc] +
        dphi[5][0]*u[ca];

    uy =
        dphi[0][1]*u[a] +
        dphi[1][1]*u[b] +
        dphi[2][1]*u[c] +
        dphi[3][1]*u[ab] +
        dphi[4][1]*u[bc] +
        dphi[5][1]*u[ca];
}

void assemble_rhs_P2(const Mesh &mesh,
                     TArray<double> &F)
{
    size_t vtx_count = mesh.vertex_count();
    size_t tri_count = mesh.triangle_count();
    const TArray<uint32_t> &idx = mesh.indices;

    // Edge numbering
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(mesh, edge_id, edge_count);

    // Inicializar RHS a cero
    for (size_t i = 0; i < F.size; ++i)
        F[i] = 0.0;

    // 6 puntos de cuadratura de orden 4 (Dunavant)
    const int Q = 6;
    double lambda[Q][3] = {
        {0.816847572980459, 0.091576213509771, 0.091576213509771},
        {0.091576213509771, 0.816847572980459, 0.091576213509771},
        {0.091576213509771, 0.091576213509771, 0.816847572980459},
        {0.108103018168070, 0.445948490915965, 0.445948490915965},
        {0.445948490915965, 0.108103018168070, 0.445948490915965},
        {0.445948490915965, 0.445948490915965, 0.108103018168070}
    };
    double w_q[Q] = {
        0.109951743655322,
        0.109951743655322,
        0.109951743655322,
        0.223381589678011,
        0.223381589678011,
        0.223381589678011
    };

    for (size_t t = 0; t < tri_count; ++t) {

        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        Vec3 A = mesh.positions[a];
        Vec3 B = mesh.positions[b];
        Vec3 C = mesh.positions[c];

        auto e0 = std::make_pair(std::min(a,b), std::max(a,b));
        auto e1 = std::make_pair(std::min(b,c), std::max(b,c));
        auto e2 = std::make_pair(std::min(c,a), std::max(c,a));

        uint32_t ab = vtx_count + edge_id[e0];
        uint32_t bc = vtx_count + edge_id[e1];
        uint32_t ca = vtx_count + edge_id[e2];

        double area = 0.5 * fabs(
            (B.x - A.x)*(C.y - A.y) -
            (C.x - A.x)*(B.y - A.y)
        );

        for (int q = 0; q < Q; ++q) {

            double l1 = lambda[q][0];
            double l2 = lambda[q][1];
            double l3 = lambda[q][2];

            double x = l1*A.x + l2*B.x + l3*C.x;
            double y = l1*A.y + l2*B.y + l3*C.y;

            double fq = exact_f(x,y);

            double phi[6];
            phi[0] = l1*(2*l1 - 1);
            phi[1] = l2*(2*l2 - 1);
            phi[2] = l3*(2*l3 - 1);
            phi[3] = 4*l1*l2;
            phi[4] = 4*l2*l3;
            phi[5] = 4*l3*l1;

            double w = w_q[q] * area;

            F[a]  += fq * phi[0] * w;
            F[b]  += fq * phi[1] * w;
            F[c]  += fq * phi[2] * w;
            F[ab] += fq * phi[3] * w;
            F[bc] += fq * phi[4] * w;
            F[ca] += fq * phi[5] * w;
        }
    }

}




void compute_errors_P1(const Mesh &mesh,
                       const TArray<double> &u_h,
                       double &L2, double &H1)
{
    L2 = 0.0;
    H1 = 0.0;

    for (size_t t = 0; t < mesh.triangle_count(); ++t) {

        int a = mesh.indices[3*t + 0];
        int b = mesh.indices[3*t + 1];
        int c = mesh.indices[3*t + 2];

        Vec3 A = mesh.positions[a];
        Vec3 B = mesh.positions[b];
        Vec3 C = mesh.positions[c];

        double area = 0.5 * fabs(
            (B.x - A.x)*(C.y - A.y) -
            (C.x - A.x)*(B.y - A.y)
        );

        // 3 puntos de cuadratura
        double lambda[3][3] = {
            {0.5, 0.5, 0.0},
            {0.5, 0.0, 0.5},
            {0.0, 0.5, 0.5}
        };

        for (int q = 0; q < 3; ++q) {

            double l1 = lambda[q][0];
            double l2 = lambda[q][1];
            double l3 = lambda[q][2];

            double x = l1*A.x + l2*B.x + l3*C.x;
            double y = l1*A.y + l2*B.y + l3*C.y;

            double u_exact = exact_u(x,y);
            double u_fem = l1*u_h[a] + l2*u_h[b] + l3*u_h[c];

            double diff = u_exact - u_fem;
            L2 += diff*diff * (area/3.0);

            // Gradiente exacto
            double ux_e, uy_e;
            exact_grad_u(x,y,ux_e,uy_e);

            // Gradiente FEM (constante en el triángulo)
            double ux_h, uy_h;
            grad_u_P1(mesh, u_h, a,b,c, ux_h, uy_h);

            double dx = ux_e - ux_h;
            double dy = uy_e - uy_h;

            H1 += (dx*dx + dy*dy) * (area/3.0);
        }
    }
}

void compute_errors_P2(const Mesh &mesh,
                       const TArray<double> &u_h,
                       double &L2, double &H1)
{
    L2 = 0.0;
    H1 = 0.0;

    size_t vtx_count = mesh.vertex_count();
    size_t tri_count = mesh.triangle_count();
    const TArray<uint32_t> &idx = mesh.indices;

    // reconstruimos edge numbering igual que en P2.cpp
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(mesh, edge_id, edge_count);

    // 6 puntos de cuadratura de orden 4 (Dunavant)
    const int Q = 6;
    double lambda[Q][3] = {
        {0.816847572980459, 0.091576213509771, 0.091576213509771},
        {0.091576213509771, 0.816847572980459, 0.091576213509771},
        {0.091576213509771, 0.091576213509771, 0.816847572980459},
        {0.108103018168070, 0.445948490915965, 0.445948490915965},
        {0.445948490915965, 0.108103018168070, 0.445948490915965},
        {0.445948490915965, 0.445948490915965, 0.108103018168070}
    };
    double w_q[Q] = {
        0.109951743655322,
        0.109951743655322,
        0.109951743655322,
        0.223381589678011,
        0.223381589678011,
        0.223381589678011
    };

    for (size_t t = 0; t < tri_count; ++t) {

        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        Vec3 A = mesh.positions[a];
        Vec3 B = mesh.positions[b];
        Vec3 C = mesh.positions[c];

        double area = 0.5 * fabs(
            (B.x - A.x)*(C.y - A.y) -
            (C.x - A.x)*(B.y - A.y)
        );

        // dofs P2: 3 vértices + 3 aristas
        auto e0 = std::make_pair(std::min(a,b), std::max(a,b));
        auto e1 = std::make_pair(std::min(b,c), std::max(b,c));
        auto e2 = std::make_pair(std::min(c,a), std::max(c,a));

        uint32_t ab = vtx_count + edge_id[e0];
        uint32_t bc = vtx_count + edge_id[e1];
        uint32_t ca = vtx_count + edge_id[e2];

        for (int q = 0; q < Q; ++q) {

            double l1 = lambda[q][0];
            double l2 = lambda[q][1];
            double l3 = lambda[q][2];

            double x = l1*A.x + l2*B.x + l3*C.x;
            double y = l1*A.y + l2*B.y + l3*C.y;

            // exacta
            double u_exact = exact_u(x,y);
            double ux_e, uy_e;
            exact_grad_u(x,y, ux_e, uy_e);

            // FEM P2
            double u_fem, ux_h, uy_h;
            eval_P2(mesh, u_h,
                    a, b, c, ab, bc, ca,
                    x, y,
                    u_fem, ux_h, uy_h);

            double du = u_exact - u_fem;
            double dx = ux_e - ux_h;
            double dy = uy_e - uy_h;

            double w = w_q[q] * area;

            L2 += du*du * w;
            H1 += (dx*dx + dy*dy) * w;
        }
    }
}



// ------------------------------
// MAIN
// ------------------------------
int main(int argc, char **argv)
{
    FEMType fem_types[2] = {FEMType::P1, FEMType::P2};
    const char* fem_names[2] = {"P1", "P2"};

    int Ns[6] = {8, 16, 32, 64, 128, 256};

    printf("N, FEM, L2_error, H1_error\n");

    for (int k = 0; k < 4; ++k) {
        int N = Ns[k];

        for (int f = 0; f < 2; ++f) {

            FEMType fem_type = fem_types[f];

            // 1. Construir malla
            Mesh mesh;
            build_unit_square_mesh(mesh, N);

            // 2. Crear solver
            PoissonSolver solver(mesh, fem_type);

            size_t vtx_count = mesh.vertex_count();

            // 3. Marcar Dirichlet en vértices (u = 0 en el borde)
            for (size_t i = 0; i < vtx_count; ++i) {
                double x = mesh.positions[i].x;
                double y = mesh.positions[i].y;

                bool boundary = (x == 0.0 || x == 1.0 || y == 0.0 || y == 1.0);

                if (boundary) {
                    solver.is_dirichlet[i] = 1;
                    solver.u_dirichlet[i] = 0.0;
                } else {
                    solver.is_dirichlet[i] = 0;
                }
            }

            // 3b. Marcar Dirichlet en dofs de arista (solo P2)
            if (fem_type == FEMType::P2) {
                std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
                uint32_t edge_count = 0;
                build_edge_numbering(mesh, edge_id, edge_count);

                for (auto &e : edge_id) {
                    uint32_t a = e.first.first;
                    uint32_t b = e.first.second;
                    uint32_t edge_dof = vtx_count + e.second;

                    double xa = mesh.positions[a].x;
                    double ya = mesh.positions[a].y;
                    double xb = mesh.positions[b].x;
                    double yb = mesh.positions[b].y;

                    bool boundary_edge =
                        (xa == 0.0 && xb == 0.0) ||
                        (xa == 1.0 && xb == 1.0) ||
                        (ya == 0.0 && yb == 0.0) ||
                        (ya == 1.0 && yb == 1.0);

                    if (boundary_edge) {
                        solver.is_dirichlet[edge_dof] = 1;
                        solver.u_dirichlet[edge_dof] = 0.0;
                    } else {
                        solver.is_dirichlet[edge_dof] = 0;
                    }
                }
            }

            if (fem_type == FEMType::P1) {

                // RHS P1 clásico: nodal sampling
                for (size_t i = 0; i < vtx_count; ++i) {
                    double x = mesh.positions[i].x;
                    double y = mesh.positions[i].y;
                    solver.f[i] = exact_f(x, y);
                }

            } else {
                // RHS P2 consistente
                assemble_rhs_P2(mesh, solver.f);
            }



            // 5. Iniciar CG (aquí dentro ya se deben aplicar las Dirichlet)
            solver.init_cg();

            // 6. Resolver Poisson
            solver.do_iterate(5000, 1e-12);

            // 7. Calcular errores
            double L2 = 0.0, H1 = 0.0;

            if (fem_type == FEMType::P1)
                compute_errors_P1(mesh, solver.u, L2, H1);
            else
                compute_errors_P2(mesh, solver.u, L2, H1);

            // 8. Imprimir fila de la tabla
            printf("%d, %s, %.10e, %.10e\n", N, fem_names[f], L2, H1);
        }
    }

    return 0;
}
