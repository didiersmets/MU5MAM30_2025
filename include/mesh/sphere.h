#include "mesh.h"

void build_face0(size_t subdiv, Mesh &mesh);
void build_face1(size_t subdiv, Mesh &mesh);
void build_face2(size_t subdiv, Mesh &mesh);
void build_face3(size_t subdiv, Mesh &mesh);
void build_face4(size_t subdiv, Mesh &mesh);
void build_face5(size_t subdiv, Mesh &mesh);

int load_sphere(Mesh &m, size_t subdiv);

/*
 * The sphere mesh is build by inscribing it into a cube, creating the cube mesh and then
 * projecting every point into the sphere.
 * The faces are ordered as follows:
 *   face:     UP DOWN FRONT RIGHT REAR LEFT
 *   face idx: 0  1    2     3     4    5c
 * The algorithm to do so works as follows:
 *  1. Every face mesh is build with the function build_facei, where the last i refers to 
       the face index.
        - The face mesh is directly added to the global mesh
        - Duplicates are not handled here
    2. Duplicates are removed
 */

/**
 * @brief Builds face UP (0) of the cube and proejcts it into the sphere
 * 
 * @param subdiv Number of vertices per edge/coordinate
 * @param mesh (Partial) global mesh
 */
template <typename T>
void build_face0<T>(size_t subdiv, Mesh &mesh) {
    T h = 1 / subdiv;  // let's hope 1/subdiv is exact floating-point...
    TVec3<T> p(-1,1,1); // top left point

    /* When building a triangle there are two possible shapes, alternating
     * Depending on which shape is in construction, index_top_point and index_prev_point
     * have slight different meanings, here shown.
     *
     * TRIANGLE 1:
     * Given the construction of a following triangle (x are vertices, A B and C are 
     * vertices of the triangle in construction:
     *     A    x    x    x  
     *     | \
     *     |  \
     *     B -- C    x    x
     * then index_top_point  - 1    =    index(A)
     *      index_prev_point        =    index(B)
     *      index_prev_point + 1    =    index(C)
     *
     * TRIANGLE 2:
     * Given the construction of a triangle (x are vertices, A B and C are vertices of 
     * the triangle in construction:
     *     x    A    x    x  
     *        / |
     *       /  |
     *     B -- C    x    x
     * then index_top_point          =    index(A)
     *      index_prev_point         =    index(B)
     *      index_prev_point + 1     =    index(C)
    */

    uint32_t index_top_point  = mesh.vertex_count() - 1;
    uint32_t index_prev_point;

    // not creating new triangles when adding the first
    for ( T x = -1; x <= 1; x += h ) {
        mesh.positions.push_back(p);
        p.x = x;
    }

    index_prev_point = mesh.vertex_count() - 1;
    // main loop (adding also triangles)
    for ( T y = 1 - h; y >= -1; y -= h ) {
        p.y = y;
        p.x = -1;

        // adding first point (no triangle creation)
        mesh.positions.push_back(p);
        // adding other points and triangles
        for ( T x = -1; x <= 1-h; x += h ) {
            // adding point
            mesh.positions.push_back(p);
            p.x = x;

            /* TRIANGLE 1 */
            mesh.indices.push_back(index_top_point-1);               // A
            mesh.indices.push_back(index_prev_point);                // B
            mesh.indices.push_back(index_prev_point + 1);   // C
            
            /* TRIANGLE 2 */
            mesh.indices.push_back(index_top_point);                 // A
            mesh.indices.push_back(index_prev_point);                // B
            mesh.indices.push_back(index_prev_point + 1);   // C

            index_top_point += 1;
            index_prev_point += 1;
        }
    }
}

void build_face1(size_t subdiv, Mesh &mesh);
void build_face2(size_t subdiv, Mesh &mesh);
void build_face3(size_t subdiv, Mesh &mesh);
void build_face4(size_t subdiv, Mesh &mesh);
void build_face5(size_t subdiv, Mesh &mesh);