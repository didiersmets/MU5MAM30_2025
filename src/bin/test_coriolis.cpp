#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "gl_utils.h"
#include "imgui/imgui.h"
#include "tiny_expr/tinyexpr.h"
#include "cube.h"
#include "logging.h"
#include "mesh.h"
#include "mesh_bounds.h"
#include "mesh_gpu.h"
#include "mesh_io.h"
#include "navier_stokes.h"
#include "ndc.h"
#include "shaders.h"
#include "sphere.h"
#include "viewer.h"

/* Viewer config */
float bgcolor[4] = {0.3, 0.3, 0.3, 1.0};
bool draw_surface = true;
bool draw_edges = false;
bool show_axes = false;  // NEW: Toggle for axis display (default off to avoid initial issues)
bool show_psi = false;  // Toggle to show psi (stream function) instead of omega
bool show_velocity_vectors = false;  // Toggle to show velocity vectors (meteorological style)
bool show_contours = false;  // Toggle to show contour lines for underlying function
bool show_omega_contours = false;  // Toggle to show omega contour lines
float axis_length = 1.5f;  // NEW: Length of axes
float velocity_scale = 0.1f;  // Scale factor for velocity arrow length
float scale_min;
float scale_max;
float mesh_deform = 0;

/* FEM interaction */
bool autoscale = true;
bool started = false;
bool one_step = false;
bool reset = false;

/* Parameters */
float lognu = -3;
float omega = 0.2f; // angular velocity
float dt = 0.0035;
double tol = 1e-6;

/* RHS expression of the PDE */
char rhs_expression[128] = "100 * x * exp(-50*x^2) * (1 + 0.5 * cos(0.05 * atan2(z, y)))";
bool rhs_show_error = false;
double rhs_x, rhs_y, rhs_z, rhs_p, rhs_t, rhs_r;
te_variable rhs_vars[] = {{"x", &rhs_x}, {"y", &rhs_y}, {"z", &rhs_z}, {"phi", &rhs_p}, {"theta", &rhs_t}, {"rand", &rhs_r}};
te_expr *te_rhs = NULL;

/* NEW: Axis rendering data */
GLuint axis_vao = 0;
GLuint axis_vbo = 0;
int axis_shader = 0;

static void syntax(char *prg_name);
static int load_mesh(Mesh &mesh, int argc, char **argv);
static void rescale_and_recenter_mesh(Mesh &mesh);
static void init_camera_for_mesh(const Mesh &mesh, Camera &camera);
static void update_all(NavierStokesSolver &solver, Mesh &mesh, GPUMesh &mesh_gpu);
static void draw_scene(const Viewer &viewer, int shader, const GPUMesh &gpu_mesh);
static void draw_gui(NavierStokesSolver &solver);
static void draw_velocity_field(const Mesh &mesh, const NavierStokesSolver &solver, const Viewer &viewer, int shader);
static void draw_contours(const Mesh &mesh, const TArray<double> &data, int shader);
static void key_cb(int key, int action, int mods, void *args);
static void get_attr_bounds(const Mesh &m, float *attr_min, float *attr_max);

/* NEW: Axis rendering functions */
static void init_axes();
static void draw_axes(const Viewer &viewer);
static void draw_axis_labels(const Viewer &viewer);
static void cleanup_axes();

void reset_solver(NavierStokesSolver &solver) {
    for (size_t i = 0; i < solver.N; ++i) {
        rhs_x = solver.m.positions[i].x;
        rhs_y = solver.m.positions[i].y;
        rhs_z = solver.m.positions[i].z;
        rhs_p = atan2(rhs_y, rhs_x);
        rhs_t = atan2(sqrt(rhs_x * rhs_x + rhs_y * rhs_y), rhs_z);
        rhs_r = (double)rand() / RAND_MAX;
        solver.omega[i] = te_eval(te_rhs);
    }
    solver.set_zero_mean(solver.omega.data);
    memset(solver.psi.data, 0, solver.N * sizeof(double));
    solver.t = 0;
}

bool new_rhs(NavierStokesSolver &solver) {
    srand((int)time(NULL));
    te_expr *test = te_compile(rhs_expression, rhs_vars, sizeof(rhs_vars) / sizeof(rhs_vars[0]), NULL);
    if (!test)
        return false;
    te_free(te_rhs);
    te_rhs = test;
    reset_solver(solver);
    return true;
}

void transfer_to_mesh(const TArray<double> &V, Mesh &m) {
    m.attr.resize(m.vertex_count());
    for (size_t i = 0; i < m.vertex_count(); ++i) {
        m.attr[i] = V[i];
    }
}

/* NEW: Initialize axis geometry */
static void init_axes() {
    // Create axis vertices (position + color)
    // Format: x, y, z, r, g, b
    float axis_vertices[] = {
        // X axis (red)
        0.0f, 0.0f, 0.0f,  1.0f, 0.0f, 0.0f,
        axis_length, 0.0f, 0.0f,  1.0f, 0.0f, 0.0f,
        // Y axis (green)
        0.0f, 0.0f, 0.0f,  0.0f, 1.0f, 0.0f,
        0.0f, axis_length, 0.0f,  0.0f, 1.0f, 0.0f,
        // Z axis (blue)
        0.0f, 0.0f, 0.0f,  0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, axis_length,  0.0f, 0.0f, 1.0f
    };

    glGenVertexArrays(1, &axis_vao);
    glGenBuffers(1, &axis_vbo);

    glBindVertexArray(axis_vao);
    glBindBuffer(GL_ARRAY_BUFFER, axis_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(axis_vertices), axis_vertices, GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Create simple shader for axes
    const char *axis_vert_src = R"(
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

    const char *axis_frag_src = R"(
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

/* Draw velocity field as meteorological wind barbs */
static void draw_velocity_field(const Mesh &mesh, const NavierStokesSolver &solver, 
                                 const Viewer &viewer, int shader)
{
    // Check if velocity has been computed (size > 0 and matches vertex count)
    if (!show_velocity_vectors || solver.velocity.size == 0 || solver.velocity.size != mesh.vertex_count()) return;
    
    // Compute velocity magnitude range for coloring
    float vel_min = 1e9, vel_max = -1e9;
    for (size_t v = 0; v < mesh.vertex_count(); ++v) {
        float vel_mag = sqrt(solver.velocity[v][0]*solver.velocity[v][0] + 
                            solver.velocity[v][1]*solver.velocity[v][1]);
        vel_min = fmin(vel_min, vel_mag);
        vel_max = fmax(vel_max, vel_mag);
    }
    if (vel_max <= vel_min) vel_max = vel_min + 1e-6;
    
    // Build vertex data for velocity arrows at vertices
    std::vector<float> vertex_data;
    std::vector<float> color_data;
    size_t vert_count = mesh.vertex_count();
    
    for (size_t v = 0; v < vert_count; ++v) {
        // Get vertex position
        Vec3f pos = mesh.positions[v];
        
        // Get velocity at this vertex
        Vec3f vel = solver.velocity[v] * velocity_scale;
        float vel_mag = sqrt(solver.velocity[v][0]*solver.velocity[v][0] + 
                            solver.velocity[v][1]*solver.velocity[v][1]);
        
        // Normalize color to [0, 1]
        float color_val = (vel_mag - vel_min) / (vel_max - vel_min);
        color_val = fmax(0.0f, fmin(1.0f, color_val));
        
        // Map to RGB: blue (low) -> green (mid) -> red (high)
        float r, g, b_col;
        if (color_val < 0.5f) {
            r = 0.0f;
            g = 2.0f * color_val;
            b_col = 1.0f - 2.0f * color_val;
        } else {
            r = 2.0f * (color_val - 0.5f);
            g = 2.0f * (1.0f - color_val);
            b_col = 0.0f;
        }
        
        // Start point (vertex position)
        vertex_data.push_back(pos.x);
        vertex_data.push_back(pos.y);
        vertex_data.push_back(pos.z);
        color_data.push_back(r);
        color_data.push_back(g);
        color_data.push_back(b_col);
        color_data.push_back(1.0f);
        
        // End point (vertex + velocity)
        vertex_data.push_back(pos.x + vel.x);
        vertex_data.push_back(pos.y + vel.y);
        vertex_data.push_back(pos.z + vel.z);
        color_data.push_back(r);
        color_data.push_back(g);
        color_data.push_back(b_col);
        color_data.push_back(1.0f);
    }
    
    if (vertex_data.empty()) return;
    
    // Create VAO and VBO for velocity vectors
    GLuint vel_vao, vel_vbo, col_vbo;
    glGenVertexArrays(1, &vel_vao);
    glGenBuffers(1, &vel_vbo);
    glGenBuffers(1, &col_vbo);
    
    glBindVertexArray(vel_vao);
    
    // Position data
    glBindBuffer(GL_ARRAY_BUFFER, vel_vbo);
    glBufferData(GL_ARRAY_BUFFER, vertex_data.size() * sizeof(float), vertex_data.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Color data
    glBindBuffer(GL_ARRAY_BUFFER, col_vbo);
    glBufferData(GL_ARRAY_BUFFER, color_data.size() * sizeof(float), color_data.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
    // Draw with shader
    glUseProgram(shader);
    glUniform1i(glGetUniformLocation(shader, "lighting"), false);
    
    const Camera &camera = viewer.camera;
    Mat4 proj = camera.view_to_clip();
    Mat4 vm = camera.world_to_view();
    glUniformMatrix4fv(glGetUniformLocation(shader, "vm"), 1, 0, &vm(0, 0));
    glUniformMatrix4fv(glGetUniformLocation(shader, "proj"), 1, 0, &proj(0, 0));
    
    glLineWidth(2.0f);
    glBindVertexArray(vel_vao);
    glDrawArrays(GL_LINES, 0, vertex_data.size() / 3);
    glBindVertexArray(0);
    glLineWidth(1.0f);
    
    // Cleanup and restore state
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDeleteBuffers(1, &vel_vbo);
    glDeleteBuffers(1, &col_vbo);
    glDeleteVertexArrays(1, &vel_vao);
    glUseProgram(0);
    while (glGetError() != GL_NO_ERROR);  // Clear any remaining errors
}

/* Draw contour lines for any function */
static void draw_contours(const Mesh &mesh, const TArray<double> &data, int shader)
{
    if (!show_contours || data.size == 0) return;
    
    // Number of contour levels
    int num_contours = 12;
    
    // Find data min/max
    double data_min = data[0];
    double data_max = data[0];
    for (size_t i = 0; i < data.size; ++i) {
        data_min = fmin(data_min, data[i]);
        data_max = fmax(data_max, data[i]);
    }
    if (data_max <= data_min) data_max = data_min + 1e-6;
    
    // Build vertex data for contour lines
    std::vector<float> vertex_data;
    std::vector<float> color_data;
    
    // For each triangle, check which contour levels it crosses
    for (size_t t = 0; t < mesh.triangle_count(); ++t) {
        uint32_t a = mesh.indices[3 * t + 0];
        uint32_t b = mesh.indices[3 * t + 1];
        uint32_t c = mesh.indices[3 * t + 2];
        
        double data_a = data[a];
        double data_b = data[b];
        double data_c = data[c];
        
        Vec3f va = mesh.positions[a];
        Vec3f vb = mesh.positions[b];
        Vec3f vc = mesh.positions[c];
        
        // For each contour level
        for (int level = 0; level < num_contours; ++level) {
            double contour_value = data_min + (level + 1.0) / (num_contours + 1.0) * (data_max - data_min);
            
            // Find edges that cross this contour
            std::vector<Vec3f> segment_points;
            
            // Check edge AB
            if ((data_a - contour_value) * (data_b - contour_value) < 0) {
                double t_param = (contour_value - data_a) / (data_b - data_a);
                Vec3f pt = va + (vb - va) * (float)t_param;
                segment_points.push_back(pt);
            }
            // Check edge BC
            if ((data_b - contour_value) * (data_c - contour_value) < 0) {
                double t_param = (contour_value - data_b) / (data_c - data_b);
                Vec3f pt = vb + (vc - vb) * (float)t_param;
                segment_points.push_back(pt);
            }
            // Check edge CA
            if ((data_c - contour_value) * (data_a - contour_value) < 0) {
                double t_param = (contour_value - data_c) / (data_a - data_c);
                Vec3f pt = vc + (va - vc) * (float)t_param;
                segment_points.push_back(pt);
            }
            
            // Draw line segment if we have exactly 2 intersection points
            if (segment_points.size() == 2) {
                vertex_data.push_back(segment_points[0].x);
                vertex_data.push_back(segment_points[0].y);
                vertex_data.push_back(segment_points[0].z);
                color_data.push_back(1.0f);  // White
                color_data.push_back(1.0f);
                color_data.push_back(1.0f);
                color_data.push_back(1.0f);
                
                vertex_data.push_back(segment_points[1].x);
                vertex_data.push_back(segment_points[1].y);
                vertex_data.push_back(segment_points[1].z);
                color_data.push_back(1.0f);  // White
                color_data.push_back(1.0f);
                color_data.push_back(1.0f);
                color_data.push_back(1.0f);
            }
        }
    }
    
    if (vertex_data.empty()) return;
    
    // Create VAO and VBO
    GLuint contour_vao, contour_vbo, contour_col_vbo;
    glGenVertexArrays(1, &contour_vao);
    glGenBuffers(1, &contour_vbo);
    glGenBuffers(1, &contour_col_vbo);
    
    glBindVertexArray(contour_vao);
    
    // Position data
    glBindBuffer(GL_ARRAY_BUFFER, contour_vbo);
    glBufferData(GL_ARRAY_BUFFER, vertex_data.size() * sizeof(float), vertex_data.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Color data
    glBindBuffer(GL_ARRAY_BUFFER, contour_col_vbo);
    glBufferData(GL_ARRAY_BUFFER, color_data.size() * sizeof(float), color_data.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    
    // Draw
    glUseProgram(shader);
    glDrawArrays(GL_LINES, 0, vertex_data.size() / 3);
    
    // Cleanup
    glDeleteBuffers(1, &contour_vbo);
    glDeleteBuffers(1, &contour_col_vbo);
    glDeleteVertexArrays(1, &contour_vao);
}

/* NEW: Draw the 3D axes */
static void draw_axes(const Viewer &viewer) {
    if (!show_axes) return;

    const Camera &camera = viewer.camera;
    
    glUseProgram(axis_shader);
    
    Mat4 proj = camera.view_to_clip();
    Mat4 vm = camera.world_to_view();
    
    GLint vm_loc = glGetUniformLocation(axis_shader, "vm");
    GLint proj_loc = glGetUniformLocation(axis_shader, "proj");
    
    if (vm_loc != -1) {
        glUniformMatrix4fv(vm_loc, 1, 0, &vm(0, 0));
    }
    if (proj_loc != -1) {
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
    while (glGetError() != GL_NO_ERROR);
}


static void draw_axis_labels(const Viewer &viewer) {
    if (!show_axes) return;
    
    const Camera &camera = viewer.camera;
    Mat4 proj = camera.view_to_clip();
    Mat4 vm = camera.world_to_view();
    Mat4 mvp = proj * vm;
    
    // Get viewport dimensions
    ImGuiIO& io = ImGui::GetIO();
    float width = io.DisplaySize.x;
    float height = io.DisplaySize.y;
    
    // Define axis endpoints in world space
    Vec3 axis_endpoints[3] = {
        Vec3(axis_length, 0, 0),  // X
        Vec3(0, axis_length, 0),  // Y
        Vec3(0, 0, axis_length)   // Z
    };
    
    const char* labels[3] = {"X", "Y", "Z"};
    ImU32 colors[3] = {
        IM_COL32(255, 0, 0, 255),    // Red for X
        IM_COL32(0, 255, 0, 255),    // Green for Y
        IM_COL32(0, 0, 255, 255)     // Blue for Z
    };
    
    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    
    for (int i = 0; i < 3; i++) {
        // Transform to clip space
        Vec3 pos = axis_endpoints[i];
        float clip_x = mvp(0,0) * pos.x + mvp(0,1) * pos.y + mvp(0,2) * pos.z + mvp(0,3);
        float clip_y = mvp(1,0) * pos.x + mvp(1,1) * pos.y + mvp(1,2) * pos.z + mvp(1,3);
        float clip_w = mvp(3,0) * pos.x + mvp(3,1) * pos.y + mvp(3,2) * pos.z + mvp(3,3);
        
        // Skip if behind camera
        if (clip_w <= 0) continue;
        
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
        draw_list->AddRectFilled(
            ImVec2(text_pos.x - 3, text_pos.y - 3),
            ImVec2(text_pos.x + text_size.x + 3, text_pos.y + text_size.y + 3),
            IM_COL32(0, 0, 0, 180)
        );
        
        // Draw text
        draw_list->AddText(text_pos, colors[i], labels[i]);
    }
}

/* NEW: Cleanup axis resources */
static void cleanup_axes() {
    if (axis_vao) {
        glDeleteVertexArrays(1, &axis_vao);
        axis_vao = 0;
    }
    if (axis_vbo) {
        glDeleteBuffers(1, &axis_vbo);
        axis_vbo = 0;
    }
    if (axis_shader) {
        glDeleteProgram(axis_shader);
        axis_shader = 0;
    }
}

int main(int argc, char **argv) {
    log_init(0);

    /* Load Mesh */
    Mesh mesh;
    if (load_mesh(mesh, argc, argv)) {
        syntax(argv[0]);
        exit(EXIT_FAILURE);
    }
    LOG_MSG("Loaded mesh.");

    rescale_and_recenter_mesh(mesh);
    LOG_MSG("Mesh rescaled and recentered.");

    /* Prepare FEM data */
    NavierStokesSolver solver(mesh);
    if (!new_rhs(solver)) {
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

    /* NEW: Initialize axes */
    init_axes();
    
    // Clear any GL errors from initialization
    while (glGetError() != GL_NO_ERROR);
    
    LOG_MSG("Axes initialized.");

    /* Main Loop */
    while (!viewer.should_close()) {
        viewer.poll_events();
        update_all(solver, mesh, gpu_mesh);

        viewer.begin_frame();
        
        // Draw scene first
        draw_scene(viewer, shader, gpu_mesh);
        
        // Draw velocity vectors if enabled
        draw_velocity_field(mesh, solver, viewer, shader);
        
        // Draw contour lines if enabled (always shows psi contours)
        if (show_contours) {
            draw_contours(mesh, solver.psi, shader);
        }
        
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

static void syntax(char *prg_name) {
    printf("Syntax : %s ($(obj_filename)| cube | sphere) [n]\n", prg_name);
    printf("  Subdivision number n must be provided in case of "
           "cube or sphere mesh.\n");
}

static int load_mesh(Mesh &mesh, int argc, char **argv) {
    int res = -1;
    if (argc > 2 && strncmp(argv[1], "cube", 4) == 0) {
        res = load_cube(mesh, atoi(argv[2]));
    } else if (argc > 2 && strncmp(argv[1], "sphere", 5) == 0) {
        res = load_sphere(mesh, atoi(argv[2]));
    } else if (argc > 1) {
        res = load_obj(argv[1], mesh);
    }
    return res;
}

static void rescale_and_recenter_mesh(Mesh &mesh) {
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

static void init_camera_for_mesh(const Mesh &mesh, Camera &camera) {
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

static void get_attr_bounds(const Mesh &m, float *attr_min, float *attr_max) {
    if (!m.vertex_count())
        return;

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

static void update_all(NavierStokesSolver &solver, Mesh &mesh, GPUMesh &gpu_mesh) {
    bool needs_upload = true;

    if (started || one_step) {
        // Compute Coriolis parameter with omega slider
        solver.time_step_coriolis(dt, pow(10, lognu), omega);
        
        // Compute velocity field from psi
        solver.compute_velocity();

        if (one_step) {
            one_step = false;
        }

        // Transfer appropriate field to mesh for visualization
        if (show_psi) {
            transfer_to_mesh(solver.psi, mesh);
        } else {
            transfer_to_mesh(solver.omega, mesh);
        }

        if (autoscale) {
            get_attr_bounds(mesh, &scale_min, &scale_max);
        }
    } else if (reset) {
        reset_solver(solver);
        if (show_psi) {
            transfer_to_mesh(solver.psi, mesh);
        } else {
            transfer_to_mesh(solver.omega, mesh);
        }
        get_attr_bounds(mesh, &scale_min, &scale_max);
        reset = false;
    } else {
        needs_upload = false;
    }

    if (needs_upload) {
        gpu_mesh.update_attr();
    }
}

static void draw_scene(const Viewer &viewer, int shader, const GPUMesh &gpu_mesh) {
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
    
    // Draw velocity field (if needed)
    // Note: draw_velocity_field accesses solver through gpu_mesh.m reference
    // We need to pass solver separately - will be called from main loop instead
}

static void draw_gui(NavierStokesSolver &solver) {
    ImGui::Begin("Controls");
    ImGui::Text("Navier Stokes solver");
    ImGui::Text("--------------------");
    ImGui::Text("Enter math expression for initial vorticity below:");
    ImGui::Text("(available variables : x, y, z, phi, theta, rand)");
    ImGui::Text("(zero mean automatically achieved by adding constant)");
    ImGui::InputText("", rhs_expression, IM_ARRAYSIZE(rhs_expression));

    if (ImGui::Button("Apply")) {
        if (!new_rhs(solver)) {
            rhs_show_error = true;
        }
        started = false;
        reset = true;
    }

    if (rhs_show_error) {
        ImGui::Begin("Error");
        ImGui::Text("Syntax error in expresion (missing * ?)");
        if (ImGui::Button("Got it!")) {
            rhs_show_error = false;
        }
        ImGui::End();
    }

    ImGui::Text(" ");
    ImGui::Text("Solution value is represented by color :");
    ImGui::Text("Red = low value, Green = mid, Blue = high.");
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

    ImGui::Text("Time : %f", solver.t);
    ImGui::Text("Scale min %.2f Scale max %.2f (Span : %g)", scale_min, scale_max, scale_max - scale_min);

    ImGui::Text(" ");
    ImGui::Text("Controls");
    ImGui::Text("--------");
    ImGui::Text("Viscosity (negative power of 10):");
    ImGui::SliderFloat("nu", &lognu, -8, 0, "10^(%.1f)");

    ImGui::Text("Angular velocity:");
    ImGui::SliderFloat("omega", &omega, 0.f, 1.f, "%.3f");

    ImGui::Text("Time step :");
    ImGui::SliderFloat("dt", &dt, 0.f, 0.01f, "%.4f");

    ImGui::Checkbox("Autoscale colors to bounds", &autoscale);
    ImGui::Checkbox("Show mesh edges", &draw_edges);
    ImGui::Checkbox("Show coordinate axes", &show_axes);
    
    ImGui::Text(" ");
    ImGui::Text("Underlying function");
    ImGui::Text("-------------------");
    if (ImGui::RadioButton("Vorticity (omega)", !show_psi)) {
        show_psi = false;
    }
    if (ImGui::RadioButton("Stream function (psi)", show_psi)) {
        show_psi = true;
    }
    
    ImGui::Text(" ");
    ImGui::Text("Overlays");
    ImGui::Text("--------");
    ImGui::Checkbox("Show contour lines", &show_contours);
    ImGui::Checkbox("Show velocity vectors", &show_velocity_vectors);
    ImGui::SliderFloat("Velocity scale", &velocity_scale, 0.01f, 0.5f);

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

static void key_cb(int key, int action, int mods, void *args) {
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

    // NEW: Toggle axes with 'A' key
    if (key == GLFW_KEY_A && action == GLFW_PRESS) {
        show_axes = !show_axes;
        return;
    }
}