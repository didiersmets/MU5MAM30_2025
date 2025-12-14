#pragma once
#include "mesh.h"
#include "duplicate_verts.h"

int load_cube(Mesh &m, size_t subdiv);
int load_cube_nested_dissect(Mesh &m, size_t subdiv);

/**
 * @brief Generates a cube mesh with overlapping vertices at edges/corners.
 *        Each face is generated separately, creating duplicate vertices.
 *        Triangles are generated with diagonal from top-left to bottom-right,
 *        all with normals pointing outward from the cube.
 * 
 * @param m The mesh to fill
 * @param subdiv Number of subdivisions per edge
 * @return int vertices per face
 */
template<typename T, bool normalize>
int load_overlapping_cube(Mesh&m, size_t subdiv);









template <typename T, bool normalize, bool invert>
void build_cube_face(Mesh& m, size_t subdiv, 
                     TVec3<T> p_corner,      // Starting corner (e.g., (-1, 1, 1) for top face)
                     TVec3<T> dir_x,         // Direction vector along X
                     TVec3<T> dir_y,         // Direction vector along Y
                     TVec3<T> normal) {      // Outward normal direction
    double h = 2.0 / subdiv;  // Step size for subdivision
    size_t start_idx = m.positions.size;  // Remember starting index for this face
    size_t vertices_per_row = subdiv + 1;
    
    // Generate all vertices for this face
    for (size_t row = 0; row <= subdiv; row++) {
        for (size_t col = 0; col <= subdiv; col++) {
            TVec3<T> p = p_corner + dir_x * (T)(col * h) + dir_y * (T)(row * h);
            if constexpr (normalize) 
                m.positions.push_back(normalized(p));
            else 
                m.positions.push_back(p);
        }
    }
    
    // Generate triangles
    for (size_t row = 0; row < subdiv; row++) {
        for (size_t col = 0; col < subdiv; col++) {
            // Current quad corners (assuming looking at the face from outside):
            // top-left, top-right, bottom-left, bottom-right
            size_t tl = start_idx + row * vertices_per_row + col;           // top-left
            size_t tr = start_idx + row * vertices_per_row + (col + 1);     // top-right
            size_t bl = start_idx + (row + 1) * vertices_per_row + col;     // bottom-left
            size_t br = start_idx + (row + 1) * vertices_per_row + (col + 1); // bottom-right
            
            // Diagonal from top-left to bottom-right
            // Triangle 1: TL, TR, BR
            // Triangle 2: TL, BR, BL
            
            // Check normal orientation: if normal points outward, use counter-clockwise
            // from the outside view. We need to check dot product with (dir_x × dir_y)
            
            // For proper winding: use right-hand rule
            // If the face normal should be dir_x × dir_y or -(dir_x × dir_y)?
            // The user wants normals pointing outward, so we need to ensure correct winding
            
            // Default: counter-clockwise when viewed from outside (normal pointing out)
            if constexpr (invert) {
                m.indices.push_back(tl);
                m.indices.push_back(br);
                m.indices.push_back(tr);
                
                m.indices.push_back(tl);
                m.indices.push_back(bl);
                m.indices.push_back(br);
            } else {
                m.indices.push_back(tl);
                m.indices.push_back(tr);
                m.indices.push_back(br);
                
                m.indices.push_back(tl);
                m.indices.push_back(br);
                m.indices.push_back(bl);
            }
        }
    }
}

template<typename T, bool normalize>
int load_overlapping_cube(Mesh&m, size_t subdiv) {
    // Cube vertices are at (±1, ±1, ±1)
    
    // Face 0: TOP (y = 1, normal = (0, 1, 0))
    // Corners: (-1, 1, -1), (1, 1, -1), (1, 1, 1), (-1, 1, 1)
    build_cube_face<T,normalize, false>(m, subdiv, 
                    TVec3<T>(-1, 1, -1),   // Starting corner
                    TVec3<T>(1, 0, 0),      // Direction X (left to right)
                    TVec3<T>(0, 0, 1),      // Direction Y (back to front)
                    TVec3<T>(0, -1, 0));     // Normal (outward)
    
    // Face 1: BOTTOM (y = -1, normal = (0, -1, 0))
    // Corners: (-1, -1, 1), (1, -1, 1), (1, -1, -1), (-1, -1, -1)
    build_cube_face<T,normalize, false>(m, subdiv,
                    TVec3<T>(-1, -1, 1),   // Starting corner
                    TVec3<T>(1, 0, 0),      // Direction X (left to right)
                    TVec3<T>(0, 0, -1),     // Direction Y (front to back)
                    TVec3<T>(0, -1, 0));    // Normal (outward)
    
    // Face 2: FRONT (z = 1, normal = (0, 0, 1))
    // Corners: (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)
    build_cube_face<T,normalize, true>(m, subdiv,
                    TVec3<T>(-1, -1, 1),   // Starting corner
                    TVec3<T>(1, 0, 0),      // Direction X (left to right)
                    TVec3<T>(0, 1, 0),      // Direction Y (down to up)
                    TVec3<T>(0, 0, 1));     // Normal (outward)
    
    // Face 3: BACK (z = -1, normal = (0, 0, -1))
    // Corners: (1, -1, -1), (-1, -1, -1), (-1, 1, -1), (1, 1, -1)
    build_cube_face<T,normalize, true>(m, subdiv,
                    TVec3<T>(1, -1, -1),   // Starting corner
                    TVec3<T>(-1, 0, 0),     // Direction X (right to left)
                    TVec3<T>(0, 1, 0),      // Direction Y (down to up)
                    TVec3<T>(0, 0, -1));    // Normal (outward)
    
    // Face 4: RIGHT (x = 1, normal = (1, 0, 0))
    // Corners: (1, -1, -1), (1, -1, 1), (1, 1, 1), (1, 1, -1)
    build_cube_face<T,normalize, false>(m, subdiv,
                    TVec3<T>(1, -1, -1),   // Starting corner
                    TVec3<T>(0, 0, 1),      // Direction X (back to front)
                    TVec3<T>(0, 1, 0),      // Direction Y (down to up)
                    TVec3<T>(1, 0, 0));     // Normal (outward)
    
    // Face 5: LEFT (x = -1, normal = (-1, 0, 0))
    // Corners: (-1, -1, 1), (-1, -1, -1), (-1, 1, -1), (-1, 1, 1)
    build_cube_face<T,normalize, false>(m, subdiv,
                    TVec3<T>(-1, -1, 1),   // Starting corner
                    TVec3<T>(0, 0, -1),     // Direction X (front to back)
                    TVec3<T>(0, 1, 0),      // Direction Y (down to up)
                    TVec3<T>(-1, 0, 0));    // Normal (outward)
    
    return subdiv + 1;
}
