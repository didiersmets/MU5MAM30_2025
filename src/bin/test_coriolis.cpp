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

/* Particle structure for velocity field visualization */
struct Particle {
    Vec3f position;           // Current position on sphere
    std::vector<Vec3f> trail; // Historical positions for fading trail
    int age;                  // Current age in time steps
    int lifespan;             // Maximum lifespan in time steps
};

typedef std::vector<Particle> ParticleList;

/* Viewer config */
float bgcolor[4] = {0.3, 0.3, 0.3, 1.0};
bool draw_surface = true;
bool draw_edges = false;
bool show_axes = false;  // NEW: Toggle for axis display (default off to avoid initial issues)
bool show_psi = false;  // Toggle to show psi (stream function) instead of omega
bool show_velocity_vectors = false;  // Toggle to show velocity vectors (meteorological style)
bool show_contours = false;  // Toggle to show contour lines for underlying function
bool show_omega_contours = false;  // Toggle to show omega contour lines
bool show_particles = false;  // Toggle to show particle tracers
float axis_length = 1.5f;  // NEW: Length of axes
float velocity_scale = 0.1f;  // Scale factor for velocity arrow length
float scale_min;
float scale_max;
float mesh_deform = 0;

/* Particle parameters */
ParticleList particles;
int particle_lifespan = 300;  // Lifespan in time steps
int particles_per_step = 2;   // Number of particles to generate per time step

/* FEM interaction */
bool autoscale = true;
bool started = false;
bool one_step = false;
bool reset = false;

/* Auto-stop */
bool enable_autostop = false;
float autostop_time = 1.0f;

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
static void draw_particles(const ParticleList &particles, const Viewer &viewer, int shader);
static void update_particles(ParticleList &particles, const NavierStokesSolver &solver, const Mesh &mesh, float dt);
static void spawn_particles(ParticleList &particles, int count, int lifespan);
static void key_cb(int key, int action, int mods, void *args);
static void get_attr_bounds(const Mesh &m, float *attr_min, float *attr_max);

/* NEW: Axis rendering functions */
static void init_axes();
static void draw_axes(const Viewer &viewer);
static void draw_axis_labels(const Viewer &viewer);
static void cleanup_axes();

/* NEW: Color bar rendering function */
static void draw_colorbar(float min_val, float max_val, const char *label);

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

/* NEW: Color mapping function matching the shader */
static void color_from_val(float u, float scale_min, float scale_max, float &r, float &g, float &b) {
    float l = (u - scale_min) / (scale_max - scale_min);
    if (l < 0.5f) {
        r = 1.0f - 2.0f * l;  // Goes from 1 to 0
        g = 2.0f * l;           // Goes from 0 to 1
        b = 0.0f;
    } else {
        l = 1.0f - l;           // Flip for second half
        r = 0.0f;
        g = 2.0f * l;           // Goes from 1 to 0
        b = 1.0f - 2.0f * l;    // Goes from 0 to 1
    }
}

/* NEW: Draw color bar function */
static void draw_colorbar(float min_val, float max_val, const char *label) {
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 220, 80), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Color Bar", nullptr);
    
    ImDrawList *draw_list = ImGui::GetWindowDrawList();
    ImVec2 canvas_pos = ImGui::GetCursorScreenPos();
    
    float bar_width = 40.0f;
    float bar_height = 350.0f;
    float bar_x = canvas_pos.x + 20.0f;
    float bar_y = canvas_pos.y + 20.0f;
    
    // Draw color gradient bar from top (max) to bottom (min)
    int num_segments = 100;
    for (int i = 0; i < num_segments; ++i) {
        // Interpolate from top (1.0) to bottom (0.0) to show max at top, min at bottom
        float t = 1.0f - (float)i / num_segments;
        
        float val = min_val + t * (max_val - min_val);
        
        float r, g, b;
        color_from_val(val, min_val, max_val, r, g, b);
        
        ImU32 col = ImGui::GetColorU32(ImVec4(r, g, b, 1.0f));
        
        float y1 = bar_y + (float)i / num_segments * bar_height;
        float y2 = bar_y + (float)(i + 1) / num_segments * bar_height;
        
        // Draw rectangle for this segment
        draw_list->AddRectFilled(
            ImVec2(bar_x, y1 ),
            ImVec2(bar_x + bar_width, y2),
            col
        );
    }
    
    // Draw border
    draw_list->AddRect(ImVec2(bar_x, bar_y), ImVec2(bar_x + bar_width, bar_y + bar_height), 
                        ImGui::GetColorU32(ImVec4(1, 1, 1, 1)), 0.0f, ImDrawCornerFlags_All, 2.0f);
    
    // Format and display title
    ImGui::Text("%s:", label);
    ImGui::Spacing();
    
    // Draw max value at top
    char max_text[32];
    snprintf(max_text, sizeof(max_text), "%.3g", max_val);
    ImGui::SetCursorScreenPos(ImVec2(bar_x + bar_width + 10, bar_y));
    ImGui::Text("%s", max_text);
    
    // Draw mid value in the middle
    char mid_text[32];
    float mid_val = (min_val + max_val) * 0.5f;
    snprintf(mid_text, sizeof(mid_text), "%.3g", mid_val);
    ImGui::SetCursorScreenPos(ImVec2(bar_x + bar_width + 10, bar_y + bar_height * 0.5f - 8.0f));
    ImGui::Text("%s", mid_text);
    
    // Draw min value at bottom
    char min_text[32];
    snprintf(min_text, sizeof(min_text), "%.3g", min_val);
    ImGui::SetCursorScreenPos(ImVec2(bar_x + bar_width + 10, bar_y + bar_height - 16.0f));
    ImGui::Text("%s", min_text);
    
    ImGui::End();
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

/* Spawn random particles on the sphere surface */
static void spawn_particles(ParticleList &particles, int count, int lifespan) {
    for (int i = 0; i < count; ++i) {
        // Generate random point on unit sphere using Fibonacci sphere algorithm
        float theta = 2.0f * M_PI * ((float)rand() / RAND_MAX);
        float phi = acos(2.0f * ((float)rand() / RAND_MAX) - 1.0f);
        
        Particle p;
        p.position.x = sin(phi) * cos(theta);
        p.position.y = sin(phi) * sin(theta);
        p.position.z = cos(phi);
        p.age = 0;
        p.lifespan = lifespan;
        p.trail.push_back(p.position);  // Initialize trail with starting position
        
        particles.push_back(p);
    }
}

/* Update particle positions using velocity field (forward Euler) */
static void update_particles(ParticleList &particles, const NavierStokesSolver &solver, const Mesh &mesh, float dt) {
    std::vector<Particle> living_particles;
    
    for (auto& particle : particles) {
        particle.age++;
        
        if (particle.age < particle.lifespan) {
            // Find velocity at particle position using barycentric interpolation
            // For now, use nearest vertex for simplicity
            float min_dist = 1e9;
            size_t nearest_vertex = 0;
            
            for (size_t v = 0; v < mesh.vertex_count(); ++v) {
                float dx = mesh.positions[v].x - particle.position.x;
                float dy = mesh.positions[v].y - particle.position.y;
                float dz = mesh.positions[v].z - particle.position.z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest_vertex = v;
                }
            }
            
            // Get velocity at nearest vertex
            Vec3f vel = solver.velocity[nearest_vertex];
            
            // Forward Euler integration
            particle.position.x += vel.x * dt;
            particle.position.y += vel.y * dt;
            particle.position.z += vel.z * dt;
            
            // Renormalize to keep particle on unit sphere
            float norm = sqrt(particle.position.x * particle.position.x + 
                             particle.position.y * particle.position.y + 
                             particle.position.z * particle.position.z);
            if (norm > 1e-6) {
                particle.position.x /= norm;
                particle.position.y /= norm;
                particle.position.z /= norm;
            }
            
            // Add current position to trail (keep trail bounded to 50 points)
            particle.trail.push_back(particle.position);
            if (particle.trail.size() > 50) {
                particle.trail.erase(particle.trail.begin());
            }
            
            living_particles.push_back(particle);
        }
    }
    
    particles = living_particles;
}

/* Draw particles as round sprites: trail spheres (fading) + small black core */
static void draw_particles(const ParticleList &particles, const Viewer &viewer, int /*shader*/) {
    if (!show_particles || particles.empty()) return;

    // Build separate buffers for trail points and particle cores
    std::vector<float> trail_pos;
    std::vector<float> trail_col;
    std::vector<float> core_pos;
    std::vector<float> core_col;

    for (const auto &p : particles) {
        // Trail samples as individual round sprites
        for (size_t i = 0; i < p.trail.size(); ++i) {
            // trail[0] is oldest, trail.back() is newest
            // age_ratio: 0.0 = oldest, 1.0 = newest
            float age_ratio = (float)i / (float)std::max<size_t>(1, p.trail.size() - 1);
            // brightness increases with recency
            float brightness = 0.2f + 0.8f * age_ratio;

            trail_pos.push_back(p.trail[i].x);
            trail_pos.push_back(p.trail[i].y);
            trail_pos.push_back(p.trail[i].z);

            // grayscale trail color with alpha fade (newer = more opaque)
            trail_col.push_back(brightness);
            trail_col.push_back(brightness);
            trail_col.push_back(brightness);
            trail_col.push_back(0.6f * age_ratio);
        }

        // Core (black) rendered on top
        core_pos.push_back(p.position.x);
        core_pos.push_back(p.position.y);
        core_pos.push_back(p.position.z);
        core_col.push_back(0.0f);
        core_col.push_back(0.0f);
        core_col.push_back(0.0f);
        core_col.push_back(1.0f);
    }

    const Camera &camera = viewer.camera;
    Mat4 proj = camera.view_to_clip();
    Mat4 vm = camera.world_to_view();

    // Simple point-sprite shader (cached)
    static GLuint point_prog = 0;
    if (!point_prog) {
        const char *vs_src = "#version 330 core\n"
            "layout(location=0) in vec3 position;\n"
            "layout(location=1) in vec4 color;\n"
            "uniform mat4 vm;\n"
            "uniform mat4 proj;\n"
            "out vec4 vColor;\n"
            "void main() {\n"
            "  gl_Position = proj * vm * vec4(position, 1.0);\n"
            "  vColor = color;\n"
            "}\n";

        const char *fs_src = "#version 330 core\n"
            "in vec4 vColor;\n"
            "out vec4 FragColor;\n"
            "void main() {\n"
            "  vec2 c = gl_PointCoord - vec2(0.5);\n"
            "  float d = length(c);\n"
            "  if (d > 0.5) discard;\n"
            "  float edge = smoothstep(0.5, 0.45, d);\n"
            "  FragColor = vec4(vColor.rgb, vColor.a * edge);\n"
            "}\n";

        GLuint vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, &vs_src, NULL);
        glCompileShader(vs);

        GLint vs_ok = GL_FALSE;
        glGetShaderiv(vs, GL_COMPILE_STATUS, &vs_ok);
        if (!vs_ok) {
            char buf[1024];
            glGetShaderInfoLog(vs, sizeof(buf), NULL, buf);
            LOG_MSG("Particle VS compile error: %s", buf);
            glDeleteShader(vs);
            // don't attempt to create program
        } else {
            GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
            glShaderSource(fs, 1, &fs_src, NULL);
            glCompileShader(fs);

            GLint fs_ok = GL_FALSE;
            glGetShaderiv(fs, GL_COMPILE_STATUS, &fs_ok);
            if (!fs_ok) {
                char buf[1024];
                glGetShaderInfoLog(fs, sizeof(buf), NULL, buf);
                LOG_MSG("Particle FS compile error: %s", buf);
                glDeleteShader(fs);
                glDeleteShader(vs);
            } else {
                point_prog = glCreateProgram();
                glAttachShader(point_prog, vs);
                glAttachShader(point_prog, fs);
                glLinkProgram(point_prog);

                GLint link_ok = GL_FALSE;
                glGetProgramiv(point_prog, GL_LINK_STATUS, &link_ok);
                if (!link_ok) {
                    char buf[1024];
                    glGetProgramInfoLog(point_prog, sizeof(buf), NULL, buf);
                    LOG_MSG("Particle shader link error: %s", buf);
                    glDeleteProgram(point_prog);
                    point_prog = 0;
                }

                glDetachShader(point_prog, vs);
                glDetachShader(point_prog, fs);
                glDeleteShader(vs);
                glDeleteShader(fs);
            }
        }
    }

    // If shader failed to compile/link, skip particle rendering
    if (!point_prog) return;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // -------- Draw trails (fading spheres) --------
    if (!trail_pos.empty()) {
        GLuint trail_vao = 0, trail_vbo = 0, trail_col_vbo = 0;
        glGenVertexArrays(1, &trail_vao);
        glGenBuffers(1, &trail_vbo);
        glGenBuffers(1, &trail_col_vbo);

        glBindVertexArray(trail_vao);

        glBindBuffer(GL_ARRAY_BUFFER, trail_vbo);
        glBufferData(GL_ARRAY_BUFFER, trail_pos.size() * sizeof(float), trail_pos.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, trail_col_vbo);
        glBufferData(GL_ARRAY_BUFFER, trail_col.size() * sizeof(float), trail_col.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glUseProgram(point_prog);
        if (glGetError() != GL_NO_ERROR) {
            LOG_MSG("glUseProgram failed for particle shader");
            glUseProgram(0);
            // cleanup generated buffers/vaos
            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);
            glDeleteBuffers(1, &trail_vbo);
            glDeleteBuffers(1, &trail_col_vbo);
            glDeleteVertexArrays(1, &trail_vao);
            glDisable(GL_BLEND);
            return;
        }
        glUniformMatrix4fv(glGetUniformLocation(point_prog, "vm"), 1, 0, &vm(0, 0));
        glUniformMatrix4fv(glGetUniformLocation(point_prog, "proj"), 1, 0, &proj(0, 0));

        // Trail points slightly larger
        glPointSize(8.0f);
        glBindVertexArray(trail_vao);
        glDrawArrays(GL_POINTS, 0, (GLsizei)(trail_pos.size() / 3));
        glBindVertexArray(0);

        glPointSize(1.0f);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDeleteBuffers(1, &trail_vbo);
        glDeleteBuffers(1, &trail_col_vbo);
        glDeleteVertexArrays(1, &trail_vao);
    }

    // -------- Draw cores (small black spheres) --------
    if (!core_pos.empty()) {
        GLuint core_vao = 0, core_vbo = 0, core_col_vbo = 0;
        glGenVertexArrays(1, &core_vao);
        glGenBuffers(1, &core_vbo);
        glGenBuffers(1, &core_col_vbo);

        glBindVertexArray(core_vao);

        glBindBuffer(GL_ARRAY_BUFFER, core_vbo);
        glBufferData(GL_ARRAY_BUFFER, core_pos.size() * sizeof(float), core_pos.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, core_col_vbo);
        glBufferData(GL_ARRAY_BUFFER, core_col.size() * sizeof(float), core_col.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glUseProgram(point_prog);
        if (glGetError() != GL_NO_ERROR) {
            LOG_MSG("glUseProgram failed for particle shader (cores)");
            glUseProgram(0);
            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);
            glDeleteBuffers(1, &core_vbo);
            glDeleteBuffers(1, &core_col_vbo);
            glDeleteVertexArrays(1, &core_vao);
            glDisable(GL_BLEND);
            return;
        }
        glUniformMatrix4fv(glGetUniformLocation(point_prog, "vm"), 1, 0, &vm(0, 0));
        glUniformMatrix4fv(glGetUniformLocation(point_prog, "proj"), 1, 0, &proj(0, 0));

        // Core points smaller so they appear as inner black dots
        glPointSize(4.0f);
        glBindVertexArray(core_vao);
        glDrawArrays(GL_POINTS, 0, (GLsizei)(core_pos.size() / 3));
        glBindVertexArray(0);
        glPointSize(1.0f);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDeleteBuffers(1, &core_vbo);
        glDeleteBuffers(1, &core_col_vbo);
        glDeleteVertexArrays(1, &core_vao);
    }

    // clear any GL errors produced here so they don't affect later draws
    while (glGetError() != GL_NO_ERROR);
    glUseProgram(0);
    glDisable(GL_BLEND);
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
        
        // Draw particles if enabled
        draw_particles(particles, viewer, shader);
        
        // Draw axes on top
        draw_axes(viewer);
        
        // Draw GUI
        draw_gui(solver);
        
        // Draw color bar with values from solver
        const char *label = show_psi ? "Psi" : "Omega";
        float bar_min = show_psi ? solver.psi_min : solver.omega_min;
        float bar_max = show_psi ? solver.psi_max : solver.omega_max;
        draw_colorbar(bar_min, bar_max, label);
        
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

/* Compute min/max values for a solver field */
static void compute_field_bounds(const TArray<double> &field, double &min_val, double &max_val) {
    if (field.size == 0) {
        min_val = 0.0;
        max_val = 0.0;
        return;
    }
    
    min_val = field[0];
    max_val = field[0];
    
    for (size_t i = 1; i < field.size; ++i) {
        if (field[i] < min_val) {
            min_val = field[i];
        } else if (field[i] > max_val) {
            max_val = field[i];
        }
    }
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

        // Auto-stop simulation if elapsed time exceeds limit
        if (enable_autostop && solver.t >= autostop_time) {
            started = false;
        }

        // Spawn particles if enabled
        if (show_particles) {
            spawn_particles(particles, particles_per_step, particle_lifespan);
            update_particles(particles, solver, mesh, dt);
        }

        // Compute min/max for both fields
        compute_field_bounds(solver.omega, solver.omega_min, solver.omega_max);
        compute_field_bounds(solver.psi, solver.psi_min, solver.psi_max);
        
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
        particles.clear();  // Clear particles on reset
        
        // Compute min/max for both fields
        compute_field_bounds(solver.omega, solver.omega_min, solver.omega_max);
        compute_field_bounds(solver.psi, solver.psi_min, solver.psi_max);
        
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
    
    ImGui::Text(" ");
    ImGui::Checkbox("Auto-stop after time", &enable_autostop);
    if (enable_autostop) {
        ImGui::InputFloat("Stop time (s)", &autostop_time, 0.1f, 0.5f, "%.2f");
    }
    
    ImGui::Text(" ");
    ImGui::Text("Time : %f", solver.t);
    ImGui::Text("Scale min %.2f Scale max %.2f (Span : %g)", scale_min, scale_max, scale_max - scale_min);

    ImGui::Text(" ");
    ImGui::Text("Controls");
    ImGui::Text("--------");
    ImGui::Text("Viscosity (negative power of 10):");
    ImGui::SliderFloat("nu", &lognu, -8, 0, "10^(%.1f)");

    ImGui::Text("Angular velocity:");
    ImGui::SliderFloat("omega", &omega, -1.f, 1.f, "%.3f");

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
    ImGui::Checkbox("Show particles", &show_particles);
    ImGui::SliderFloat("Velocity scale", &velocity_scale, 0.01f, 0.5f);
    ImGui::SliderInt("Particle lifespan", &particle_lifespan, 10, 500);
    ImGui::SliderInt("Particles per step", &particles_per_step, 1, 20);

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