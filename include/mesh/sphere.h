#include "mesh.h"
#include <vector>
#include "utils_debug.h"

#define DEBUG 1

struct FaceBoundary {
    std::vector<size_t> up;
    std::vector<size_t> down;
    std::vector<size_t> right;
    std::vector<size_t> left;
};

int load_sphere(Mesh &m, size_t subdiv);


/**
 * @brief Builds a generic cube face 2x2 of the cube and projects it into the sphere.
 * 
 * @param h Dimension of one segment
 * @param mesh (Partial) global mesh
 *
 * @param T basic data type (typically float or double)
 * @param face face number (UP=0, DOWN=1, FRONT=2, RIGHT=3, BACK=4, LEFT=5)
 */
template <typename T, int face>
void build_face(double h, Mesh &mesh, FaceBoundary& old_edges, FaceBoundary& new_edges);



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


 template <typename T>
 void _build_face1(double h, Mesh &mesh, int row, int col, TVec3<T> p, FaceBoundary& new_edges) {

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

    uint32_t index_top_point  = mesh.vertex_count();
    uint32_t index_cur_point  = mesh.vertex_count();


    /* ============================== HANDLING FIRST ROW ============================= */
    new_edges.left.push_back(index_cur_point);
    for ( T j = 1; j >= -1; j -= h ) {
        p[col] = j;
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.up.push_back(index_cur_point);

        index_cur_point++;
    }
    new_edges.right.push_back(index_cur_point-1);

    /* ============================== HANDLING CENTRAL ROWs ============================= */
    for ( T i = 1 - h; i >= -1+h; i -= h ) {
        p[row] = i;
        p[col] = 1;


        // adding first point (no triangle creation)
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.left.push_back(index_cur_point);

        index_top_point++;
        index_cur_point++;

        // adding other points and triangles
        for ( T j = 1-h; j >= -1; j -= h ) {
            // adding point
            p[col] = j;
            #ifdef DEBUG
            mesh.positions.push_back(p);
            #else
            mesh.positions.push_back(normalized(p));
            #endif

            /* TRIANGLE 1 */
            mesh.indices.push_back(index_top_point-1);               // A
            mesh.indices.push_back(index_cur_point-1);                // B
            mesh.indices.push_back(index_cur_point);            // C
            
            /* TRIANGLE 2 */
            mesh.indices.push_back(index_top_point-1);                 // A
            mesh.indices.push_back(index_cur_point);                // B
            mesh.indices.push_back(index_top_point);            // C

            index_top_point++;
            index_cur_point++;
        }
        new_edges.right.push_back(index_cur_point - 1);
    }

    /* ============================== HANDLING FINAL ROW ============================= */
    p[row] = -1;
    p[col] = 1;


    // adding first point (no triangle creation)
    #ifdef DEBUG
    mesh.positions.push_back(p);
    #else
    mesh.positions.push_back(normalized(p));
    #endif
    new_edges.left.push_back(index_cur_point);
    new_edges.down.push_back(index_cur_point);

    index_top_point++;
    index_cur_point++;

    // adding other points and triangles
    for ( T j = 1-h; j >= -1; j -= h ) {
        // adding point
        p[col] = j;
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.down.push_back(index_cur_point);

        /* TRIANGLE 1 */
        mesh.indices.push_back(index_top_point-1);               // A
        mesh.indices.push_back(index_cur_point-1);                // B
        mesh.indices.push_back(index_cur_point);            // C
        
        /* TRIANGLE 2 */
        mesh.indices.push_back(index_top_point-1);                 // A
        mesh.indices.push_back(index_cur_point);                // B
        mesh.indices.push_back(index_top_point);            // C

        index_top_point++;
        index_cur_point++;
    }
    new_edges.right.push_back(index_cur_point - 1);
 }

 template <typename T>
 void _build_face0(double h, Mesh &mesh, int row, int col, TVec3<T> p, FaceBoundary& new_edges) {


    /* When building a triangle there are two possible shapes, alternating
     * Depending on which shape is in construction, index_top_point and index_prev_point
     * have slight different meanings, here shown.
     *
     * TRIANGLE 1:
     * Given the construction of a following triangle (x are vertices, A B and C are 
     * vertices of the triangle in construction:
     *     A    B    x    x  
     *     | \
     *     |  \
     *     C -- D
     * then index_top_point  - 1    =    index(A)
     *      index_prev_point        =    index(C)
     *      index_prev_point + 1    =    index(D)
     *
     * TRIANGLE 2:
     * Given the construction of a triangle (x are vertices, A B and C are vertices of 
     * the triangle in construction:
     *     A -- B    x    x  
     *       \  |
     *        \ |
     *     C -- D    x    x
     * then index_top_point  - 1     =    index(A)
     *      index_prev_point + 1     =    index(D)
     *      index_top_point          =    index(B)
    */

    uint32_t index_top_point  = mesh.vertex_count();
    uint32_t index_cur_point = mesh.vertex_count();

    /* ============================== HANDLING FIRST ROW ============================= */
    new_edges.left.push_back(index_cur_point);
    for ( T j = -1; j <= 1; j += h ) {
        p[col] = j;
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.up.push_back(index_cur_point);

        index_cur_point++;
    }
    new_edges.right.push_back(index_cur_point-1);

    /* ============================== HANDLING CENTRAL ROWS ============================= */
    for ( T i = 1 - h; i >= -1+h; i -= h ) {
        p[row] = i;
        p[col] = -1;


        // adding first point (no triangle creation)
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.left.push_back(index_cur_point);

        index_top_point++;
        index_cur_point++;

        // adding other points and triangles
        for ( T j = -1+h; j <= 1; j += h ) {
            // adding point
            p[col] = j;
            #ifdef DEBUG
            mesh.positions.push_back(p);
            #else
            mesh.positions.push_back(normalized(p));
            #endif

            /* TRIANGLE 1 */
            mesh.indices.push_back(index_top_point - 1);               // A
            mesh.indices.push_back(index_cur_point - 1);                // C
            mesh.indices.push_back(index_cur_point);            // D
            
            /* TRIANGLE 2 */
            mesh.indices.push_back(index_top_point - 1);             // A
            mesh.indices.push_back(index_cur_point);            // D
            mesh.indices.push_back(index_top_point);                 // B

            index_top_point++;
            index_cur_point++;

        }
        new_edges.right.push_back(index_cur_point - 1);
    }

    /* ============================== HANDLING LAST ROWS ============================= */
    p[row] = -1;
    p[col] = -1;

    // adding first point (no triangle creation)
    #ifdef DEBUG
    mesh.positions.push_back(p);
    #else
    mesh.positions.push_back(normalized(p));
    #endif
    new_edges.left.push_back(index_cur_point);
    new_edges.down.push_back(index_cur_point);

    index_top_point++;
    index_cur_point++;

    // adding other points and triangles
    for ( T j = -1+h; j <= 1; j += h ) {
        // adding point
        p[col] = j;
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.down.push_back(index_cur_point);

        /* TRIANGLE 1 */
        mesh.indices.push_back(index_top_point - 1);               // A
        mesh.indices.push_back(index_cur_point - 1);                // C
        mesh.indices.push_back(index_cur_point);            // D
        
        /* TRIANGLE 2 */
        mesh.indices.push_back(index_top_point - 1);             // A
        mesh.indices.push_back(index_cur_point);            // D
        mesh.indices.push_back(index_top_point);                 // B

        index_top_point++;
        index_cur_point++;

    }
    new_edges.right.push_back(index_cur_point - 1);

 }

 template <typename T>
 void _build_face2(double h, Mesh &mesh, int row, int col, TVec3<T> p, FaceBoundary& old_edges, FaceBoundary& new_edges) {

    uint32_t index_top_point = mesh.vertex_count();
    uint32_t index_cur_point = mesh.vertex_count(); // points to the first place free

    /* ============================== HANDLING FIRST 2 ROWS ============================= */
    p[row] = 1 - h;
    p[col] = -1;
    #ifdef DEBUG
    mesh.positions.push_back(p);
    #else
    mesh.positions.push_back(normalized(p));
    #endif
    new_edges.left.push_back(old_edges.up[0]);
    new_edges.left.push_back(index_cur_point);

    index_cur_point++;

    for ( int j = 1; j <= 2.0/h; j++ ) {
        p[col] = -1 + j*h;
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        
        /* TRIANGLE 1 */
        mesh.indices.push_back(old_edges.up[j-1]);               // A
        mesh.indices.push_back(index_cur_point - 1);                // C
        mesh.indices.push_back(index_cur_point);            // D
        
        /* TRIANGLE 2 */
        mesh.indices.push_back(old_edges.up[j-1]);             // A
        mesh.indices.push_back(index_cur_point);            // D
        mesh.indices.push_back(old_edges.up[j]);                 // B

        index_cur_point++;
    }

    new_edges.right.push_back(old_edges.up[2.0/h]);
    new_edges.right.push_back(index_cur_point-1);

    /* ============================== HANDLING CENTRAL ROWS ============================= */
    for ( T i = 1 - 2* h; i >= -1 + h; i -= h ) {
        p[row] = i;
        p[col] = -1;


        // adding first point (no triangle creation)
        #ifdef DEBUG
        mesh.positions.push_back(p);
        #else
        mesh.positions.push_back(normalized(p));
        #endif
        new_edges.left.push_back(index_cur_point);

        index_top_point++;
        index_cur_point++;

        // adding other points and triangles
        for ( T j = -1+h; j <= 1; j += h ) {
            // adding point
            p[col] = j;
            #ifdef DEBUG
            mesh.positions.push_back(p);
            #else
            mesh.positions.push_back(normalized(p));
            #endif

            /* TRIANGLE 1 */
            mesh.indices.push_back(index_top_point - 1);               // A
            mesh.indices.push_back(index_cur_point - 1);                // C
            mesh.indices.push_back(index_cur_point);            // D
            
            /* TRIANGLE 2 */
            mesh.indices.push_back(index_top_point - 1);             // A
            mesh.indices.push_back(index_cur_point + 1);            // D
            mesh.indices.push_back(index_top_point);                 // B

            index_top_point++;
            index_cur_point++;

        }

        new_edges.right.push_back(index_cur_point-1);
    }

    /* ============================== HANDLING FINAL ROW ============================= */
    new_edges.left.push_back(old_edges.down[0]);

    index_top_point++;
    index_cur_point++;
    for ( int j = 1; j <= 2.0/h; j++ ) {
        /* TRIANGLE 1 */
        mesh.indices.push_back(index_top_point - 1);               // A
        mesh.indices.push_back(old_edges.down[j-1]);                // C
        mesh.indices.push_back(old_edges.down[j]);            // D
        
        /* TRIANGLE 2 */
        mesh.indices.push_back(index_top_point - 1);             // A
        mesh.indices.push_back(old_edges.down[j]);            // D
        mesh.indices.push_back(index_top_point);                 // B

        index_top_point++;
    }
    new_edges.right.push_back(old_edges.down[2.0/h]);

 }

/**
 * @brief Builds a generic cube face 2x2 of the cube and projects it into the sphere.
 * 
 * @param h Dimension of one segment
 * @param mesh (Partial) global mesh
 *
 * @param T basic data type (typically float or double)
 * @param face face number (UP=0, DOWN=1, FRONT=2, RIGHT=3, BACK=4, LEFT=5)
 */
template <typename T, int face>
void build_face(double h, Mesh &mesh, FaceBoundary& old_edges, FaceBoundary& new_edges) {
    int row, col;
    //TVec3<T> p(-1,1,1); // top left point
    TVec3<T>p(0,0,0);
    if constexpr ( face == 0 ) {     // UP
        row = 1; col = 0;  // y, x
        p.x = -1; p.y = 1; p.z = 1;
        _build_face0<T>(h, mesh, row, col, p, new_edges);
    } else if    ( face == 1 ) {     // DOWN
        row = 1; col = 0;  // y, x
        p.x = 1; p.y = 1; p.z = -1;
        _build_face1<T>(h, mesh, row, col, p, new_edges);
    } else if    ( face == 2 ) {     // FRONT
        row = 2; col = 0;  // y, x
        p.x = -1; p.y = -1; p.z = 1;
        _build_face2<T>(h, mesh, row, col, p, old_edges, new_edges);
    } else if    ( face == 4 ) {     // BACK
        row = 2; col = 0;  // y, x
        p.x = 1; p.y = 1; p.z = 1;
        // _build_face4<T>(h, mesh, row, col, p, old_edges, new_edges);
    } else if    ( face == 3 ) {     // RIGHT
        row = 2; col = 1;  // y, x
        p.x = 1; p.y = -1; p.z = 1;
        // _build_face3<T>(h, mesh, row, col, p, old_edges, new_edges);
    } else if    ( face == 5 ) {     // LEFT
        row = 2; col = 1;  // y, x
        p.x = -1; p.y = 1; p.z = -1;
        // _build_face5<T>(h, mesh, row, col, p, old_edges, new_edges);
    }

}


int load_sphere(Mesh &m, size_t subdiv) {
    // ======================== MESH CREATION =========================
    double h = 2.0 / subdiv;
    FaceBoundary in, out_up, out_down, out_front, out_right, out_back, out_left;

    build_face<float,0>(h, m, in, out_up);
    print_vector(out_up.up);
    print_vector(out_up.down);
    print_vector(out_up.right);
    print_vector(out_up.left);
    printf("\n");
    build_face<float,1>(h, m, in, out_down);
    print_vector(out_down.up);
    print_vector(out_down.down);
    print_vector(out_down.right);
    print_vector(out_down.left);
    printf("\n");

    in.up   = out_up.down;
    in.down = out_down.down;
    in.left.clear();
    in.right.clear();
    build_face<float,2>(h, m, in, out_front);
    print_vector(out_front.up);
    print_vector(out_front.down);
    print_vector(out_front.right);
    print_vector(out_front.left);

    // in.up   = out_up.right;
    // in.down = out_down.left;
    // in.left = out_front.right;
    // in.right.clear();
    // build_face<float,3>(subdiv, m, in, out_right);

    // in.up   = out_up.up; 
    // in.down = out_down.up;
    // in.left = out_right.right;
    // build_face<float,4>(subdiv, m, in, out_back);

    // in.up   = out_up.left; 
    // in.down = out_down.right;
    // in.left = out_back.right;
    // in.right= out_front.left;
    // build_face<float,5>(subdiv, m, in, out_left);
    
    return 0;
}
