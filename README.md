Github repo for the students of the course MU5MAM30 : Implementation of the Finite Element Method.
Tutorial classes on Wednesday 1.30pm to 4.30pm in 16-26 401

To get write permission access to this repo, first clone it using the _token_ that 
will be given to you in the first lecture. All of your work shall be uploaded into
a separate branch named according to your family name. In seqence :

1. Clone repository : \
*git clone https://REPLACE_WITH_TOKEN_HERE@github.com/didiersmets/MU5MAM30_2025.git*
2. Move into it : \
*cd MU5MAM30_2025*
3. Create and switch to a new branch named after your family name : \
*git checkout -B YOUR_FAMILY_NAME*


![test_viewer](./data/test_viewer.jpg?raw=true "Test FEM viewer")

---

# MU5MAM30 — Minimal Graph Surface Solver

This project implements a finite-element solver for the **minimal graph problem** on 2D domains.
Given Dirichlet boundary data `f` on `∂Ω`, the solver finds the surface `u : Ω → ℝ` that minimises the area functional

```
E(u) = ∫_Ω √(1 + |∇u|²) dx
```

Two nonlinear iterative schemes are provided: **Newton** and **Picard** (referred to as *Picardi* in the code).

---

## Repository structure

```
.
├── CMakeLists.txt              # Build configuration
├── configure.sh                # CMake configuration helper
├── make.sh                     # Build helper
├── shaders/                    # GLSL shaders for the viewer
├── extern/                     # Third-party libraries (glfw, imgui, tiny_expr)
├── include/
│   ├── common/                 # Math primitives (Vec2, Vec3, Mat3, Mat4, …)
│   ├── fem/                    # FEM headers: P1.h, minimal_graph.h, …
│   ├── linalg/                 # Linear solvers (conjugate_gradient.h)
│   ├── matrix/                 # Matrix types (sparse_matrix.h, fem_matrix.h)
│   └── mesh/                   # Mesh types and generators (mesh.h, square.h, disk.h)
└── src/
    ├── bin/
    │   ├── test_minimal_surface.cpp           # Interactive viewer (main entry point)
    │   └── test_minimal_surface_experiments.cpp  # Batch experiment runner
    ├── fem/
    │   ├── P1.cpp              # P1 FEM assembly (stiffness, mass, weighted variants)
    │   └── minimal_graph.cpp   # Newton and Picard solvers (MinimalGraphSolver)
    ├── linalg/
    │   └── conjugate_gradient.cpp
    ├── matrix/
    │   ├── sparse_matrix.cpp   # CSR matrix and operator()
    │   └── fem_matrix.cpp      # FEMatrix matrix-vector products
    └── mesh/
        ├── square.cpp          # Structured square mesh generator
        ├── disk.cpp            # Disk mesh generator
        └── adjacency.cpp       # Vertex-triangle adjacency
```

---

## How to build

### 1. Configure (once)

```bash
bash configure.sh
```

This runs CMake and creates two build directories:

| Directory       | Build type |
|-----------------|------------|
| `build/debug`   | Debug      |
| `build/release` | Release    |

### 2. Compile

```bash
bash make.sh
```

Or build a single target manually:

```bash
make --no-print-directory -C build/debug test_minimal_surface -j4
```

The two main executables produced are:

| Executable                           | Description                       |
|--------------------------------------|-----------------------------------|
| `test_minimal_surface`               | Interactive OpenGL viewer         |
| `test_minimal_surface_experiments`   | Headless batch experiment runner  |

---

## How to launch `test_minimal_surface`

```
./test_minimal_surface (square | disk) <n> [size1] [size2] [newton|picardi]
```

| Argument          | Meaning                                                      |
|-------------------|--------------------------------------------------------------|
| `square` / `disk` | Domain shape                                                 |
| `n`               | Subdivision level (cells per side / radial rings)            |
| `size1`           | Half-side length (square) or radius (disk). Default: `1.0`   |
| `size2`           | Second half-side length (square only). Default: `1.0`        |
| `newton`          | Use Newton iteration (default)                               |
| `picardi`         | Use Picard iteration                                         |

**Examples:**

```bash
# Square domain [0,1]² with 20 subdivisions, Newton solver
./test_minimal_surface square 20 1.0 1.0 newton

# Disk of radius 2 with 30 radial rings, Picard solver
./test_minimal_surface disk 30 2.0 picardi
```

### Interactive controls

| Control                | Action                    |
|------------------------|---------------------------|
| GUI → **Start / Stop** | Run / pause the solver    |
| GUI → **One step**     | Advance one iteration     |
| GUI → **Reset**        | Reset solution to BCs     |
| `S`                    | Toggle surface rendering  |
| `E`                    | Toggle edge rendering     |
| Mouse drag             | Orbit camera              |
| `Ctrl` + drag          | Zoom                      |
| `Shift` + drag         | Translate                 |

---

## Main function — `src/bin/test_minimal_surface.cpp`

`main()` performs the following steps:

1. **Parse arguments** — `parse_solver_mode` selects Newton or Picard; `load_mesh` calls the appropriate mesh generator.
2. **Rescale and recenter** — `rescale_and_recenter_mesh` normalises coordinates to `[-1, 1]³` for display.
3. **Construct the solver** — `MinimalGraphSolver solver(mesh, test_f)` where `test_f` is the Dirichlet boundary function. The current default is the Scherk-type exact minimal surface:
   ```
   f(x,y) = (1/α) · ln( cos(αx) / cos(αy) )    with α = π/4
   ```
4. **Initialise viewer and GPU resources** — OpenGL context, GLSL shaders, GPU mesh upload.
5. **Enter the render loop** — on each frame, `update_all()` calls `do_iterate_Newton()` or `do_iterate_Picardi()` for `iter_per_frame` steps, copies `solver.u` into mesh vertex positions (`z` coordinate), and re-uploads to the GPU.

---

## Mesh initialisation

### Square mesh — [`src/mesh/square.cpp`](src/mesh/square.cpp)

```cpp
void build_square_mesh(Mesh *m, size_t N, double a, double b);
```

Builds a structured `(N+1) × (N+1)` vertex grid over `[0,a] × [0,b]`, split into `2N²` triangles
(each grid square is divided along its diagonal).
All four sides are detected and stored in `m->boundary`.

### Disk mesh — [`src/mesh/disk.cpp`](src/mesh/disk.cpp)

```cpp
void build_disk_mesh(Mesh *m, size_t N, double R);
```

Builds a disk of radius `R` with `N` radial rings.
The outermost ring of vertices is stored in `m->boundary`.

---

## Linear system assembly

### Common data structures

| Symbol | Type | Description |
|--------|------|-------------|
| `P` | `CSRPattern` | Sparsity pattern of the stiffness matrix (upper triangular, symmetric). Built once by `build_P1_CSRPattern`. |
| `S` | `CSRMatrix` | Assembled stiffness matrix, sharing `row_start` / `col` pointers with `P`. |
| `q` | `TArray<double>` | Per-triangle denominator `qₜ = 1 / √(1 + |∇uₜ|²)`, updated every iteration. |

---

### Step 1 — Compute the denominator `q`

```cpp
double MinimalGraphSolver::compute_denominator(TArray<double> &den,
                                               const TArray<double> &u);
// src/fem/minimal_graph.cpp
```

For each triangle `T` with vertices `(a, b, c)`:

1. Computes local stiffness entries `S_loc[0..5]` from the edge vectors `AB`, `AC`.
2. Evaluates `|∇u|²_T = uᵀ S_loc u / area_T`.
3. Stores `den[T] = 1 / √(1 + |∇u|²_T)`.

Returns the current surface area `∑_T area_T / den[T]` (used as the energy in the Newton line search).

---

### Newton iteration — [`src/fem/minimal_graph.cpp`](src/fem/minimal_graph.cpp) — `do_iterate_Newton`

At each Newton step the linear system `J(uₖ) δu = −F(uₖ)` is solved for `δu`, then
`u_{k+1} = u_k + α δu` with `α` chosen by Armijo backtracking.

#### Step 2 — Jacobian assembly (left-hand side)

```cpp
void build_P1_stiffness_matrix_NS(const Mesh &m, const CSRPattern &P,
                                  CSRMatrix &S, const double *den,
                                  const double *u, double area);
// src/fem/P1.cpp
```

For each triangle `T`, the local Jacobian block is:

```
J_loc[i,j] = qₜ · S_loc[i,j]  −  qₜ³ · (∇u · ∇φᵢ)(∇u · ∇φⱼ) / area_T
```

where `S_loc` is the standard P1 stiffness matrix for `T` and `qₜ = den[T]`.
Boundary rows are penalised: `J(bᵢ, bᵢ) = 1e30`.

#### Step 3 — Residual / right-hand side

```cpp
void build_P1_rhs_NS(const Mesh &m, const double *den, const double *u,
                     TArray<double> &rhs);
// src/fem/P1.cpp
```

Computes `F(u)ᵢ = ∑_T qₜ · (S_loc · u)ᵢ` (the discrete gradient of the area functional).
Boundary entries are set to zero so the CG increment `δu` is zero at boundary nodes.

---

### Picard iteration — [`src/fem/minimal_graph.cpp`](src/fem/minimal_graph.cpp) — `do_iterate_Picardi`

The Picard scheme freezes `q` at the previous iterate and solves the symmetric linear system:

```
A(uₖ) uₖ₊₁ = rhs
```

Dirichlet boundary conditions are enforced by the penalty method:
`A(bᵢ, bᵢ) = 1e30` and `rhs[bᵢ] = f(bᵢ) × 1e30`.

#### Step 2 — Weighted stiffness matrix (left-hand side)

```cpp
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P,
                               CSRMatrix &S, const double *den);
// src/fem/P1.cpp
```

Assembles `A(u)ᵢⱼ = ∑_T den[T] · S_loc[i,j]` — the standard P1 stiffness matrix with each
triangle's local contribution scaled by `qₜ = den[T]`.

#### Step 3 — Right-hand side

The right-hand side is zero for all interior nodes.
For boundary nodes, `b[bᵢ] = f(bᵢ)` is set in `clear_solution` and then multiplied by
`1e30` at the start of `do_iterate_Picardi`, so that the penalty system enforces `u[bᵢ] ≈ f(bᵢ)`.

---

## Key source files at a glance

| File | Responsibility |
|------|----------------|
| [`src/fem/minimal_graph.cpp`](src/fem/minimal_graph.cpp) | `MinimalGraphSolver`: Newton + Picard drivers, denominator computation, solver state |
| [`src/fem/P1.cpp`](src/fem/P1.cpp) | CSR stiffness / mass assembly, Newton Jacobian, Picard weighted stiffness, RHS |
| [`src/mesh/square.cpp`](src/mesh/square.cpp) | Structured square mesh + boundary extraction |
| [`src/mesh/disk.cpp`](src/mesh/disk.cpp) | Disk mesh + boundary extraction |
| [`src/linalg/conjugate_gradient.cpp`](src/linalg/conjugate_gradient.cpp) | Preconditioned and plain CG solver |
| [`src/matrix/sparse_matrix.cpp`](src/matrix/sparse_matrix.cpp) | CSR `operator()`, matrix-vector product |
| [`include/fem/minimal_graph.h`](include/fem/minimal_graph.h) | `MinimalGraphSolver` struct declaration |
| [`include/fem/P1.h`](include/fem/P1.h) | Declarations for all P1 assembly functions |
