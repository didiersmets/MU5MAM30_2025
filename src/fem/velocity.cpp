#include "navier_stokes.h"
#include "vec3.h"
#include <vector>
#include <cmath>

void NavierStokesSolver::compute_velocity()
{
  // Compute velocity field as the perpendicular gradient of psi at VERTICES
  // u = ∇⊥ psi = normal × ∇psi (on the tangent plane)
  // Since psi is P1 Lagrange, compute by averaging gradient contributions from adjacent triangles
  
  size_t vert_count = m.vertex_count();
  velocity.resize(vert_count);
  
  // Initialize velocity to zero
  for (size_t v = 0; v < vert_count; ++v) {
    velocity[v] = Vec3::Zero;
  }
  
  // Accumulate gradient contributions from each triangle to its vertices
  std::vector<double> vertex_weights(vert_count, 0.0);
  
  size_t tri_count = m.triangle_count();
  for (size_t t = 0; t < tri_count; ++t) {
    // Get vertex indices
    uint32_t a = m.indices[3 * t + 0];
    uint32_t b = m.indices[3 * t + 1];
    uint32_t c = m.indices[3 * t + 2];
    
    // Get vertex positions
    Vec3f pa = m.positions[a];
    Vec3f pb = m.positions[b];
    Vec3f pc = m.positions[c];
    
    // Convert to double for computation
    Vec3d A = {(double)pa[0], (double)pa[1], (double)pa[2]};
    Vec3d B = {(double)pb[0], (double)pb[1], (double)pb[2]};
    Vec3d C = {(double)pc[0], (double)pc[1], (double)pc[2]};
    
    // Compute triangle area
    Vec3d AB = B - A;
    Vec3d AC = C - A;
    Vec3d cross_prod = cross(AB, AC);
    double twice_area = norm(cross_prod);
    
    if (twice_area < 1e-14) {
      continue;  // Skip degenerate triangles
    }
    
    // Triangle normal (surface normal for gradient computation)
    Vec3d tri_normal = cross_prod * (1.0 / twice_area);
    
    // For P1 Lagrange elements, gradient is constant on the element.
    double factor = 1.0 / twice_area;
    
    // Edge vectors in 3D
    Vec3d BC = C - B;  // Opposite to vertex a
    Vec3d CA = A - C;  // Opposite to vertex b
    Vec3d AB_edge = B - A;  // Opposite to vertex c
    
    // Project edge vectors onto tangent plane at triangle
    auto project_to_tangent = [&](Vec3d e) -> Vec3d {
      double e_dot_n = e[0] * tri_normal[0] + e[1] * tri_normal[1] + e[2] * tri_normal[2];
      return e - tri_normal * e_dot_n;
    };
    
    Vec3d BC_tan = project_to_tangent(BC);
    Vec3d CA_tan = project_to_tangent(CA);
    Vec3d AB_tan = project_to_tangent(AB_edge);
    
    // For P1 basis: ∇φ_i is perpendicular to opposite edge (in tangent plane)
    // Perpendicular in 3D: v⊥ = normal × v
    auto perpendicular_tangent = [&](Vec3d e_tan) -> Vec3d {
      return cross(tri_normal, e_tan);
    };
    
    Vec3d grad_phi_a = perpendicular_tangent(BC_tan) * factor;
    Vec3d grad_phi_b = perpendicular_tangent(CA_tan) * factor;
    Vec3d grad_phi_c = perpendicular_tangent(AB_tan) * factor;
    
    // Weight by triangle area for averaging
    double area_weight = twice_area * 0.5;
    
    // Accumulate gradient contribution to each vertex
    Vec3d grad_psi = grad_phi_a * psi[a] + grad_phi_b * psi[b] + grad_phi_c * psi[c];
    Vec3d vel_contrib = cross(tri_normal, grad_psi);
    
    // Cast to Vec3 (float) for accumulation
    Vec3 vel_contrib_f = Vec3{(float)vel_contrib[0], (float)vel_contrib[1], (float)vel_contrib[2]};
    float weight_f = (float)(area_weight / 3.0);
    
    velocity[a] = velocity[a] + vel_contrib_f * weight_f;
    velocity[b] = velocity[b] + vel_contrib_f * weight_f;
    velocity[c] = velocity[c] + vel_contrib_f * weight_f;
    
    vertex_weights[a] += area_weight / 3.0;
    vertex_weights[b] += area_weight / 3.0;
    vertex_weights[c] += area_weight / 3.0;
  }
  
  // Normalize by accumulated weight
  for (size_t v = 0; v < vert_count; ++v) {
    if (vertex_weights[v] > 1e-14) {
      velocity[v] = velocity[v] * (float)(1.0 / vertex_weights[v]);
    } else {
      velocity[v] = Vec3::Zero;
    }
  }
}
