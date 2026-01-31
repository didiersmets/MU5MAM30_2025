#include "common/gl_utils.h"
#include "common/logging.h"
#include "fem/navier_stokes.h"
#include "imgui/imgui.h"
#include "mesh/cube.h"
#include "mesh/mesh.h"
#include "mesh/mesh_bounds.h"
#include "mesh/mesh_gpu.h"
#include "mesh/mesh_io.h"
#include "mesh/sphere.h"
#include "tiny_expr/tinyexpr.h"
#include "viewer/ndc.h"
#include "viewer/shaders.h"
#include "viewer/viewer.h"

#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Viewer config */
float bgcolor[4]   = {0.3, 0.3, 0.3, 1.0};
bool  draw_surface = true;
bool  draw_edges   = false;
bool  show_axes    = false;  // NEW: Toggle for axis display (default off to avoid initial issues)
float axis_length  = 1.5f;   // NEW: Length of axes
float scale_min;
float scale_max;
float mesh_deform = 0;

/* FEM interaction */
bool autoscale = true;
bool started   = false;
bool one_step  = false;
bool reset     = false;

/* Parameters */
float  lognu       = -3.7;
float  dt          = 0.002;
float  denominator = 1;  // denominator for angular velocity
double tol         = 1e-6;

/* RHS expression of the PDE */
char        rhs_expression[128] = "100 * x * exp(-50*x^2) * (1 + 0.5 * cos(0.05 * atan2(z, y)))";
bool        rhs_show_error      = false;
double      rhs_x, rhs_y, rhs_z, rhs_p, rhs_t, rhs_r;
te_variable rhs_vars[] = {{"x", &rhs_x},
                          {"y", &rhs_y},
                          {"z", &rhs_z},
                          {"phi", &rhs_p},
                          {"theta", &rhs_t},
                          {"rand", &rhs_r}};
te_expr*    te_rhs     = NULL;

/* NEW: Axis rendering data */
GLuint axis_vao    = 0;
GLuint axis_vbo    = 0;
int    axis_shader = 0;

static void syntax(char* prg_name);
static int  load_mesh(Mesh& mesh, int argc, char** argv);
static void rescale_and_recenter_mesh(Mesh& mesh);
static void init_camera_for_mesh(const Mesh& mesh, Camera& camera);
static void update_all(NavierStokesSolver& solver, Mesh& mesh, GPUMesh& mesh_gpu);
static void draw_scene(const Viewer& viewer, int shader, const GPUMesh& gpu_mesh);
static void draw_gui(NavierStokesSolver& solver);
static void key_cb(int key, int action, int mods, void* args);
static void get_attr_bounds(const Mesh& m, float* attr_min, float* attr_max);

/* NEW: Axis rendering functions */
static void init_axes();
static void draw_axes(const Viewer& viewer);
static void draw_axis_labels(const Viewer& viewer);
static void cleanup_axes();

void reset_solver(NavierStokesSolver& solver)
{
  for (size_t i = 0; i < solver.N; ++i)
  {
    rhs_x           = solver.m.positions[i].x;
    rhs_y           = solver.m.positions[i].y;
    rhs_z           = solver.m.positions[i].z;
    rhs_p           = atan2(rhs_y, rhs_x);
    rhs_t           = atan2(sqrt(rhs_x * rhs_x + rhs_y * rhs_y), rhs_z);
    rhs_r           = (double) rand() / RAND_MAX;
    solver.omega[i] = te_eval(te_rhs);
  }
  solver.set_zero_mean(solver.omega.data);
  memset(solver.psi.data, 0, solver.N * sizeof(double));
  solver.t = 0;
}

bool new_rhs(NavierStokesSolver& solver)
{
  srand((int) time(NULL));
  te_expr* test =
    te_compile(rhs_expression, rhs_vars, sizeof(rhs_vars) / sizeof(rhs_vars[0]), NULL);
  if (!test)
    return false;
  te_free(te_rhs);
  te_rhs = test;
  reset_solver(solver);
  return true;
}

void transfer_to_mesh(const TArray<double>& V, Mesh& m)
{
  m.attr.resize(m.vertex_count());
  for (size_t i = 0; i < m.vertex_count(); ++i)
  {
    m.attr[i] = V[i];
  }
}

/* NEW: Initialize axis geometry */
static void init_axes()
{
  // Create axis vertices (position + color)
  // Format: x, y, z, r, g, b
  float axis_vertices[] = {// X axis (red)
                           0.0f,
                           0.0f,
                           0.0f,
                           1.0f,
                           0.0f,
                           0.0f,
                           axis_length,
                           0.0f,
                           0.0f,
                           1.0f,
                           0.0f,
                           0.0f,
                           // Y axis (green)
                           0.0f,
                           0.0f,
                           0.0f,
                           0.0f,
                           1.0f,
                           0.0f,
                           0.0f,
                           axis_length,
                           0.0f,
                           0.0f,
                           1.0f,
                           0.0f,
                           // Z axis (blue)
                           0.0f,
                           0.0f,
                           0.0f,
                           0.0f,
                           0.0f,
                           1.0f,
                           0.0f,
                           0.0f,
                           axis_length,
                           0.0f,
                           0.0f,
                           1.0f};

  glGenVertexArrays(1, &axis_vao);
  glGenBuffers(1, &axis_vbo);

  glBindVertexArray(axis_vao);
  glBindBuffer(GL_ARRAY_BUFFER, axis_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(axis_vertices), axis_vertices, GL_STATIC_DRAW);

  // Position attribute
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*) 0);
  glEnableVertexAttribArray(0);

  // Color attribute
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*) (3 * sizeof(float)));
  glEnableVertexAttribArray(1);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  // Create simple shader for axes
  const char* axis_vert_src = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aColor;
        
        uniform mat4 vm;
        uniform mat4 proj;
        
        out vec3 Color;
        
        void main() {
            gl_Position = proj * vm * vec4(aPos, 1.0);
            Color = aColor;
        }
    )";

  const char* axis_frag_src = R"(
        #version 330 core
        in vec3 Color;
        out vec4 FragColor;
        
        void main() {
            FragColor = vec4(Color, 1.0);
        }
    )";

  // Compile shaders
  GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vert_shader, 1, &axis_vert_src, NULL);
  glCompileShader(vert_shader);

  GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(frag_shader, 1, &axis_frag_src, NULL);
  glCompileShader(frag_shader);

  axis_shader = glCreateProgram();
  glAttachShader(axis_shader, vert_shader);
  glAttachShader(axis_shader, frag_shader);
  glLinkProgram(axis_shader);

  glDeleteShader(vert_shader);
  glDeleteShader(frag_shader);
}

/* NEW: Draw the 3D axes */
static void draw_axes(const Viewer& viewer)
{
  if (!show_axes)
    return;

  const Camera& camera = viewer.camera;

  glUseProgram(axis_shader);

  Mat4 proj = camera.view_to_clip();
  Mat4 vm   = camera.world_to_view();

  GLint vm_loc   = glGetUniformLocation(axis_shader, "vm");
  GLint proj_loc = glGetUniformLocation(axis_shader, "proj");

  if (vm_loc != -1)
  {
    glUniformMatrix4fv(vm_loc, 1, 0, &vm(0, 0));
  }
  if (proj_loc != -1)
  {
    glUniformMatrix4fv(proj_loc, 1, 0, &proj(0, 0));
  }

  glBindVertexArray(axis_vao);

  // Disable line smoothing to ensure thick lines work
  glDisable(GL_LINE_SMOOTH);

  // Set thick line width (hardcoded to 8.0)
  glLineWidth(8.0f);

  glDrawArrays(GL_LINES, 0, 6);

  // Reset line width
  glLineWidth(1.0f);
  glBindVertexArray(0);

  // Clear any errors that may have occurred during axis rendering
  while (glGetError() != GL_NO_ERROR)
    ;
}

/* NEW: Draw axis labels using ImGui */
static void draw_axis_labels(const Viewer& viewer)
{
  if (!show_axes)
    return;

  const Camera& camera = viewer.camera;
  Mat4          proj   = camera.view_to_clip();
  Mat4          vm     = camera.world_to_view();
  Mat4          mvp    = proj * vm;

  // Get viewport dimensions
  ImGuiIO& io     = ImGui::GetIO();
  float    width  = io.DisplaySize.x;
  float    height = io.DisplaySize.y;

  // Define axis endpoints in world space
  Vec3 axis_endpoints[3] = {
    Vec3(axis_length, 0, 0),  // X
    Vec3(0, axis_length, 0),  // Y
    Vec3(0, 0, axis_length)   // Z
  };

  const char* labels[3] = {"X", "Y", "Z"};
  ImU32       colors[3] = {
    IM_COL32(255, 0, 0, 255),  // Red for X
    IM_COL32(0, 255, 0, 255),  // Green for Y
    IM_COL32(0, 0, 255, 255)   // Blue for Z
  };

  ImDrawList* draw_list = ImGui::GetForegroundDrawList();

  for (int i = 0; i < 3; i++)
  {
    // Transform to clip space
    Vec3  pos    = axis_endpoints[i];
    float clip_x = mvp(0, 0) * pos.x + mvp(0, 1) * pos.y + mvp(0, 2) * pos.z + mvp(0, 3);
    float clip_y = mvp(1, 0) * pos.x + mvp(1, 1) * pos.y + mvp(1, 2) * pos.z + mvp(1, 3);
    float clip_w = mvp(3, 0) * pos.x + mvp(3, 1) * pos.y + mvp(3, 2) * pos.z + mvp(3, 3);

    // Skip if behind camera
    if (clip_w <= 0)
      continue;

    // Convert to NDC
    float ndc_x = clip_x / clip_w;
    float ndc_y = clip_y / clip_w;

    // Convert to screen space
    float screen_x = (ndc_x * 0.5f + 0.5f) * width;
    float screen_y = (1.0f - (ndc_y * 0.5f + 0.5f)) * height;

    // Draw label with background for better visibility
    ImVec2 text_pos(screen_x + 10, screen_y - 10);
    ImVec2 text_size = ImGui::CalcTextSize(labels[i]);

    // Draw background rectangle
    draw_list->AddRectFilled(ImVec2(text_pos.x - 3, text_pos.y - 3),
                             ImVec2(text_pos.x + text_size.x + 3, text_pos.y + text_size.y + 3),
                             IM_COL32(0, 0, 0, 180));

    // Draw text
    draw_list->AddText(text_pos, colors[i], labels[i]);
  }
}

/* NEW: Cleanup axis resources */
static void cleanup_axes()
{
  if (axis_vao)
  {
    glDeleteVertexArrays(1, &axis_vao);
    axis_vao = 0;
  }
  if (axis_vbo)
  {
    glDeleteBuffers(1, &axis_vbo);
    axis_vbo = 0;
  }
  if (axis_shader)
  {
    glDeleteProgram(axis_shader);
    axis_shader = 0;
  }
}

int main(int argc, char** argv)
{
  log_init(0);

  /* Load Mesh */
  Mesh mesh;
  if (load_mesh(mesh, argc, argv))
  {
    syntax(argv[0]);
    exit(EXIT_FAILURE);
  }
  LOG_MSG("Loaded mesh.");

  rescale_and_recenter_mesh(mesh);
  LOG_MSG("Mesh rescaled and recentered.");

  /* Prepare FEM data */
  NavierStokesSolver solver(mesh);
  if (!new_rhs(solver))
  {
    LOG_MSG("Error loading rhs (expression flawed ?).");
    exit(EXIT_FAILURE);
  }
  transfer_to_mesh(solver.omega, mesh);
  get_attr_bounds(mesh, &scale_min, &scale_max);
  LOG_MSG("Prepared FEM data.");

  /* Get an OpenGL context through a viewer app. */
  Viewer viewer;
  init_camera_for_mesh(mesh, viewer.camera);
  viewer.init("Navier Stokes 2D solver (vorticity formulation)");
  viewer.register_key_callback({key_cb, NULL});
  viewer.mouse.set_double_click_time(-1);
  LOG_MSG("Viewer initialized.");

  /* Prepare GPU data */
  const char* vert_shader = "./shaders/fem.vert";
  const char* frag_shader = "./shaders/fem.frag";
  int         shader      = create_shader(vert_shader, frag_shader);
  if (!shader)
  {
    exit(EXIT_FAILURE);
  }
  LOG_MSG("Shader initialized.");

  GPUMesh gpu_mesh;
  gpu_mesh.m = &mesh;
  gpu_mesh.upload();

  /* NEW: Initialize axes */
  init_axes();

  // Clear any GL errors from initialization
  while (glGetError() != GL_NO_ERROR)
    ;

  LOG_MSG("Axes initialized.");

  /* Main Loop */
  while (!viewer.should_close())
  {
    viewer.poll_events();
    update_all(solver, mesh, gpu_mesh);

    viewer.begin_frame();

    // Draw scene first
    draw_scene(viewer, shader, gpu_mesh);

    // Draw axes on top
    draw_axes(viewer);

    // Draw GUI
    draw_gui(solver);

    // Draw axis labels (after GUI so they appear on top)
    draw_axis_labels(viewer);

    viewer.end_frame();
  }

  /* NEW: Cleanup axes */
  cleanup_axes();

  viewer.fini();
  log_fini();

  return (EXIT_SUCCESS);
}

static void syntax(char* prg_name)
{
  printf("Syntax : %s ($(obj_filename)| cube | sphere) [n]\n", prg_name);
  printf("  Subdivision number n must be provided in case of "
         "cube or sphere mesh.\n");
}

static int load_mesh(Mesh& mesh, int argc, char** argv)
{
  int res = -1;
  if (argc > 2 && strncmp(argv[1], "cube", 4) == 0)
  {
    res = load_cube(mesh, atoi(argv[2]));
  }
  else if (argc > 2 && strncmp(argv[1], "sphere", 5) == 0)
  {
    res = load_sphere(mesh, atoi(argv[2]));
  }
  else if (argc > 1)
  {
    res = load_obj(argv[1], mesh);
  }
  return res;
}

static void rescale_and_recenter_mesh(Mesh& mesh)
{
  Aabb  bbox         = compute_mesh_bounds(mesh);
  Vec3  model_center = (bbox.min + bbox.max) * 0.5f;
  Vec3  model_extent = (bbox.max - bbox.min);
  float model_size   = max(model_extent);

  if (model_size == 0)
  {
    printf("Warning : Mesh is empty or reduced to a point.\n");
    model_size = 1;
  }

  for (size_t i = 0; i < mesh.vertex_count(); ++i)
  {
    mesh.positions[i] -= model_center;
    mesh.positions[i] /= (model_size / 2);
  }
}

static void init_camera_for_mesh(const Mesh& mesh, Camera& camera)
{
  Aabb  bbox         = compute_mesh_bounds(mesh);
  Vec3  model_center = (bbox.min + bbox.max) * 0.5f;
  Vec3  model_extent = (bbox.max - bbox.min);
  float model_size   = max(model_extent);

  if (model_size == 0)
  {
    printf("Warning : Mesh is empty or reduced to a point.\n");
    model_size = 1;
  }

  camera.set_target(model_center);
  Vec3 start_pos = (model_center + 2.f * Vec3(0, 0, model_size));
  camera.set_position(start_pos);
  camera.set_near(0.01 * model_size);
  camera.set_far(100 * model_size);
}

static void get_attr_bounds(const Mesh& m, float* attr_min, float* attr_max)
{
  if (!m.vertex_count())
    return;

  float min = m.attr[0];
  float max = min;

  for (size_t i = 1; i < m.vertex_count(); ++i)
  {
    if (m.attr[i] < min)
    {
      min = m.attr[i];
    }
    else if (m.attr[i] > max)
    {
      max = m.attr[i];
    }
  }

  *attr_min = min;
  *attr_max = max;
}

static void update_all(NavierStokesSolver& solver, Mesh& mesh, GPUMesh& gpu_mesh)
{
  bool needs_upload = true;

  if (started || one_step)
  {
    // Calculate omega as 10^-denominator.
    // If denominator = 4.13, we obtain the actual Earth angular velocity (~7.29e-5 rad/s).
    float omega_value = 1.0f / pow(10.0f, denominator);

    solver.time_step_coriolis(dt, pow(10, lognu), (double) omega_value);

    if (one_step)
    {
      one_step = false;
    }

    transfer_to_mesh(solver.omega, mesh);

    if (autoscale)
    {
      get_attr_bounds(mesh, &scale_min, &scale_max);
    }
  }
  else if (reset)
  {
    reset_solver(solver);
    transfer_to_mesh(solver.omega, mesh);
    get_attr_bounds(mesh, &scale_min, &scale_max);
    reset = false;
  }
  else
  {
    needs_upload = false;
  }

  if (needs_upload)
  {
    gpu_mesh.update_attr();
  }
}

static void draw_scene(const Viewer& viewer, int shader, const GPUMesh& gpu_mesh)
{
  glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], bgcolor[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  const Camera& camera = viewer.camera;

  glUseProgram(shader);

  Mat4 proj       = camera.view_to_clip();
  Mat4 vm         = camera.world_to_view();
  Vec3 camera_pos = camera.get_position();

  glUniformMatrix4fv(glGetUniformLocation(shader, "vm"), 1, 0, &vm(0, 0));
  glUniformMatrix4fv(glGetUniformLocation(shader, "proj"), 1, 0, &proj(0, 0));
  glUniform3fv(glGetUniformLocation(shader, "camera_pos"), 1, &camera_pos[0]);
  glUniform1f(glGetUniformLocation(shader, "scale_min"), scale_min);
  glUniform1f(glGetUniformLocation(shader, "scale_max"), scale_max);
  glUniform1f(glGetUniformLocation(shader, "deform"), mesh_deform);

  if (draw_surface)
  {
    glEnable(GL_POLYGON_OFFSET_FILL);
    float offset = reversed_z ? -1.f : 1.f;
    glPolygonOffset(offset, offset);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glUniform1i(glGetUniformLocation(shader, "lighting"), true);
    gpu_mesh.draw();
  }

  if (draw_edges)
  {
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(0.f, 0.f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glUniform1i(glGetUniformLocation(shader, "lighting"), false);
    gpu_mesh.draw();
  }
}

static void draw_gui(NavierStokesSolver& solver)
{
  ImGui::Begin("Controls");
  ImGui::Text("Navier Stokes solver");
  ImGui::Text("--------------------");
  ImGui::Text("Enter math expression for initial vorticity below:");
  ImGui::Text("(available variables : x, y, z, phi, theta, rand)");
  ImGui::Text("(zero mean automatically achieved by adding constant)");
  ImGui::InputText("", rhs_expression, IM_ARRAYSIZE(rhs_expression));

  if (ImGui::Button("Apply"))
  {
    if (!new_rhs(solver))
    {
      rhs_show_error = true;
    }
    started = false;
    reset   = true;
  }

  if (rhs_show_error)
  {
    ImGui::Begin("Error");
    ImGui::Text("Syntax error in expresion (missing * ?)");
    if (ImGui::Button("Got it!"))
    {
      rhs_show_error = false;
    }
    ImGui::End();
  }

  ImGui::Text(" ");
  ImGui::Text("Solution value is represented by color :");
  ImGui::Text("Red = low value, Green = mid, Blue = high.");
  ImGui::Text(" ");

  if (ImGui::Button("Start"))
  {
    started = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("Stop"))
  {
    started = false;
  }
  ImGui::SameLine();
  if (ImGui::Button("One step"))
  {
    if (!started)
    {
      one_step = true;
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Reset"))
  {
    reset = true;
  }

  ImGui::Text("Time : %f", solver.t);
  ImGui::Text("Scale min %.2f Scale max %.2f (Span : %g)",
              scale_min,
              scale_max,
              scale_max - scale_min);

  ImGui::Text(" ");
  ImGui::Text("Controls");
  ImGui::Text("--------");
  ImGui::Text("Viscosity (negative power of 10):");
  ImGui::SliderFloat("nu", &lognu, -8, 0, "10^(%.1f)");

  // --- OMEGA SLIDER SETUP ---
  // A range from 0 to 6 provides optimal control:
  // 0    = Fast rotation (1 rad/s)
  // 4.13 = Real Earth rotation (~7.29e-5 rad/s)
  // 6    = Nearly imperceptible rotation
  ImGui::Text("Angular velocity (10^-n rad/s):");
  ImGui::SliderFloat("omega exponent", &denominator, 0.0f, 6.0f, "n = %.1f");

  // Display the actual numerical value for user convenience
  float current_omega = 1.0f / std::pow(10.0f, denominator);
  ImGui::Text("Current Omega: %.2e rad/s", (double) current_omega);

  ImGui::Text("Time step :");
  ImGui::SliderFloat("dt", &dt, 0.f, 0.01f, "%.4f");

  ImGui::Checkbox("Autoscale colors to bounds", &autoscale);
  ImGui::Checkbox("Show mesh edges", &draw_edges);
  ImGui::Checkbox("Show coordinate axes", &show_axes);  // NEW: Toggle for axes

  ImGui::Text("Artificially deform mesh according to omega :");
  ImGui::Text("(may help visualize oscillations)");
  ImGui::SliderFloat(" ", &mesh_deform, 0.f, 1.f);

  ImGui::Text(" ");
  ImGui::Text("Number of DOF : %zu", solver.N);
  float fps = ImGui::GetIO().Framerate;
  ImGui::Text("Average framerate : %.1f FPS", fps);

  ImGui::Text(" ");
  ImGui::Text("Mouse :");
  ImGui::Text("Click + drag : orbit");
  ImGui::Text("Click + CTRL + drag : zoom in/out");
  ImGui::Text("Click + SHIFT + drag : translate");

  ImGui::Text(" ");
  ImGui::Text("Axes: Red=X, Green=Y, Blue=Z");  // NEW: Legend

  ImGui::End();
}

static void key_cb(int key, int action, int mods, void* args)
{
  (void) mods;
  (void) args;

  if (key == GLFW_KEY_S && action == GLFW_PRESS)
  {
    draw_surface = !draw_surface;
    return;
  }

  if (key == GLFW_KEY_E && action == GLFW_PRESS)
  {
    draw_edges = !draw_edges;
    return;
  }

  // NEW: Toggle axes with 'A' key
  if (key == GLFW_KEY_A && action == GLFW_PRESS)
  {
    show_axes = !show_axes;
    return;
  }
}