#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>

#include "gl_utils.h"
#include "imgui/imgui.h"
#include "logging.h"
#include "cube.h"
#include "mesh.h"
#include "square.h"
#include "disk.h"
#include "minimal_graph.h"
#include "mesh_gpu.h"
#include "mesh_bounds.h"
#include "ndc.h"
#include "shaders.h"
#include "sphere.h"
#include "viewer.h"

/* Viewer config */
float bgcolor[4] = {0.3f, 0.3f, 0.3f, 1.0f};
bool draw_surface = true;
bool draw_edges = false;
float scale_min;
float scale_max;
float mesh_deform = 0.0f;

/* Solver interaction */
bool autoscale = true;
bool started = false;
bool one_step = false;
bool reset = false;
int iter_per_frame = 1;
enum class SolverMode {
    Newton,
    Picardi
};
SolverMode solver_mode = SolverMode::Newton;

/*
static double test_f(const Vec2d &pos)
{
    double x = pos.x;
    double y = pos.y;
    return 1.0 / 3.0 * (pow(x - 2.0, 3) + 3.0 * (x - 2.0) * pow(y - 2.0, 2));
}
    */

/*
static double test_f(const Vec2d &pos)
{
    double x = pos.x;
    double y = pos.y;
    return 1.0 / 3.0 * (x * x * x - 3.0 * x * y * y);
}
*/

static double test_f(const Vec2d &pos)
{
    constexpr double alpha = M_PI / 4.0;
    double x = pos.x;
    double y = pos.y;
    return (1.0 / alpha) * std::log(std::cos(alpha * x) / std::cos(alpha * y));
}


static void syntax(char *prg_name)
{
	printf("Syntax : %s (square | disk) [n] [size1] [size2] [newton|picardi]\n", prg_name);
	printf("         Disk uses size1 as radius; square uses size1 and size2.\n");
	printf("         The last optional argument selects the solver mode.\n");
}

static SolverMode parse_solver_mode(int argc, char **argv)
{
    const char *mode = argv[argc - 1];
    if (strcmp(mode, "newton") == 0) {
        return SolverMode::Newton;
    }
    if (strcmp(mode, "picardi") == 0 || strcmp(mode, "picard") == 0) {
        return SolverMode::Picardi;
    }
    return SolverMode::Newton;
}

static int load_mesh(Mesh &mesh, int argc, char **argv)
{
    if (argc < 3) {
        return -1;
    }
    const size_t subdiv = static_cast<size_t>(atoi(argv[2]));
    if (strcmp(argv[1], "disk") == 0) {
        double radius = (argc > 3) ? atof(argv[3]) : 1.0;
        if (radius <= 0.0) {
            return -1;
        }
        build_disk_mesh(&mesh, subdiv, radius);
        return 0;
    }
    if (strcmp(argv[1], "square") == 0) {
        double a = (argc > 3) ? atof(argv[3]) : 1.0;
        double b = (argc > 4) ? atof(argv[4]) : 1.0;
        if (a <= 0.0 || b <= 0.0) {
            return -1;
        }
        build_square_mesh(&mesh, subdiv, a, b);
        return 0;
    }
    return -1;
}

static void transfer_to_mesh(const TArray<double> &V, Mesh &m)
{
    m.attr.resize(m.vertex_count());
    for (size_t i = 0; i < m.vertex_count(); ++i) {
        m.attr[i] = static_cast<float>(V[i]);
        m.positions[i].z = static_cast<float>(V[i]);
    }
}

static void rescale_and_recenter_mesh(Mesh &mesh)
{
    Aabb bbox = compute_mesh_bounds(mesh);
    Vec3 model_center = (bbox.min + bbox.max) * 0.5f;
    Vec3 model_extent = (bbox.max - bbox.min);
    float model_size = max(model_extent);
    if (model_size == 0) {
        printf("Warning : Mesh is empty or reduced to a point.\n");
        model_size = 1;
    }
    for (size_t i = 0; i < mesh.vertex_count(); ++i) {
        mesh.positions[i] -= model_center;
        mesh.positions[i] /= (model_size / 2);
    }
}

static void init_camera_for_mesh(const Mesh &mesh, Camera &camera)
{
    Aabb bbox = compute_mesh_bounds(mesh);
    Vec3 model_center = (bbox.min + bbox.max) * 0.5f;
    Vec3 model_extent = (bbox.max - bbox.min);
    float model_size = max(model_extent);
    if (model_size == 0) {
        printf("Warning : Mesh is empty or reduced to a point.\n");
        model_size = 1;
    }
    camera.set_target(model_center);
    Vec3 start_pos = (model_center + 2.f * Vec3(0, 0, model_size));
    camera.set_position(start_pos);
    camera.set_near(0.01 * model_size);
    camera.set_far(100 * model_size);
}

static void get_attr_bounds(const Mesh &m, float *attr_min, float *attr_max)
{
    if (!m.vertex_count()) {
        return;
    }
    float min = m.attr[0];
    float max = min;
    for (size_t i = 1; i < m.vertex_count(); ++i) {
        if (m.attr[i] < min) {
            min = m.attr[i];
        } else if (m.attr[i] > max) {
            max = m.attr[i];
        }
    }
    *attr_min = min;
    *attr_max = max;
}

static void draw_scene(const Viewer &viewer, int shader, const GPUMesh &gpu_mesh)
{
    glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], bgcolor[3]);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const Camera &camera = viewer.camera;
    glUseProgram(shader);
    Mat4 proj = camera.view_to_clip();
    Mat4 vm = camera.world_to_view();
    Vec3 camera_pos = camera.get_position();
    glUniformMatrix4fv(glGetUniformLocation(shader, "vm"), 1, 0, &vm(0, 0));
    glUniformMatrix4fv(glGetUniformLocation(shader, "proj"), 1, 0, &proj(0, 0));
    glUniform3fv(glGetUniformLocation(shader, "camera_pos"), 1, &camera_pos[0]);
    glUniform1f(glGetUniformLocation(shader, "scale_min"), scale_min);
    glUniform1f(glGetUniformLocation(shader, "scale_max"), scale_max);
    glUniform1f(glGetUniformLocation(shader, "deform"), mesh_deform);

    if (draw_surface) {
        glEnable(GL_POLYGON_OFFSET_FILL);
        float offset = reversed_z ? -1.f : 1.f;
        glPolygonOffset(offset, offset);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glUniform1i(glGetUniformLocation(shader, "lighting"), true);
        gpu_mesh.draw();
    }
    if (draw_edges) {
        glDisable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(0.f, 0.f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glUniform1i(glGetUniformLocation(shader, "lighting"), false);
        gpu_mesh.draw();
    }
}

static void update_all(MinimalGraphSolver &solver, Mesh &mesh, GPUMesh &gpu_mesh)
{
    bool needs_upload = true;
    if (started || one_step) {
		if (solver_mode == SolverMode::Newton) {
			solver.do_iterate_Newton(iter_per_frame, 1e-12, 0.1);
		} else {
			solver.do_iterate_Picardi(iter_per_frame, 1e-12);
		}
        transfer_to_mesh(solver.u, mesh);
        if (autoscale) {
            get_attr_bounds(mesh, &scale_min, &scale_max);
        }
        if (one_step) {
            one_step = false;
        }
    } else if (reset) {
		solver.clear_solution(solver_mode == SolverMode::Newton);
        transfer_to_mesh(solver.u, mesh);
        get_attr_bounds(mesh, &scale_min, &scale_max);
        reset = false;
    } else {
        needs_upload = false;
    }
    if (needs_upload) {
        gpu_mesh.update_positions();
        gpu_mesh.update_attr();
    }
    if (solver.converged) {
        started = false;
    }
}

static void draw_gui(MinimalGraphSolver &solver)
{
    ImGui::Begin("Controls");
    ImGui::Text("Minimal graph solver");
    ImGui::Text("----------------------");
    ImGui::Text("Solver : %s", solver_mode == SolverMode::Newton ? "Newton" : "Picardi");
    ImGui::Text("Iterate : %zu", solver_mode == SolverMode::Newton ? solver.iterate_N : solver.iterate_P);
    ImGui::Text("Converged : %s", solver.converged ? "yes" : "no");
    ImGui::Text("Number of DOF : %zu", solver.N);
    ImGui::Text("Scale min %.2f Scale max %.2f", scale_min, scale_max);

    ImGui::Text(" ");
    if (ImGui::Button("Start")) {
        started = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Stop")) {
        started = false;
    }
    ImGui::SameLine();
    if (ImGui::Button("One step")) {
        if (!started) {
            one_step = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
        reset = true;
    }

    ImGui::Text(" ");
    ImGui::Checkbox("Autoscale", &autoscale);
    ImGui::Checkbox("Show edges", &draw_edges);
    ImGui::Text("Iterations per frame :");
    ImGui::DragInt(" ", &iter_per_frame, 1, 1, 20);
    ImGui::Text(" ");
    ImGui::Text("Mouse :");
    ImGui::Text("Click + drag : orbit");
    ImGui::Text("Click + CTRL + drag : zoom in/out");
    ImGui::Text("Click + SHIFT + drag : translate");
    ImGui::End();
}

static void key_cb(int key, int action, int mods, void *args)
{
    (void)mods;
    (void)args;
    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        draw_surface = !draw_surface;
        return;
    }
    if (key == GLFW_KEY_E && action == GLFW_PRESS) {
        draw_edges = !draw_edges;
        return;
    }
}

int main(int argc, char **argv)
{
    log_init(0);

    solver_mode = parse_solver_mode(argc, argv);

    Mesh mesh;
    if (load_mesh(mesh, argc, argv)) {
        syntax(argv[0]);
        exit(EXIT_FAILURE);
    }
    LOG_MSG("Loaded mesh.");
    rescale_and_recenter_mesh(mesh);
    LOG_MSG("Mesh rescaled and recentered.");

    MinimalGraphSolver solver(mesh, test_f);
	solver.clear_solution(solver_mode == SolverMode::Newton);
	transfer_to_mesh(solver.u, mesh);
    get_attr_bounds(mesh, &scale_min, &scale_max);

    Viewer viewer;
    init_camera_for_mesh(mesh, viewer.camera);
    viewer.init("Minimal graph solver (Mesh init)");
    viewer.register_key_callback({key_cb, NULL});
	LOG_MSG("Viewer initialized.");

    const char *vert_shader = "./shaders/fem.vert";
    const char *frag_shader = "./shaders/fem.frag";
    int shader = create_shader(vert_shader, frag_shader);
    if (!shader) {
		exit(EXIT_FAILURE);
	}
	LOG_MSG("Shader initialized.");

    GPUMesh gpu_mesh;
    gpu_mesh.m = &mesh;
    gpu_mesh.upload();

    while (!viewer.should_close()) {
        viewer.poll_events();
        update_all(solver, mesh, gpu_mesh);
        viewer.begin_frame();
        draw_scene(viewer, shader, gpu_mesh);
        draw_gui(solver);
        viewer.end_frame();
    }

    viewer.fini();
    log_fini();
    return EXIT_SUCCESS;
}
