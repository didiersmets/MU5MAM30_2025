#include "mesh/mesh.h"
#include <map>

struct _face_mesh {
  TArray<Vec3> positions;
  TArray<uint32_t> indices;
  TArray<int> vertices_of_type[3];
};


void _load_diamond(_face_mesh &m, size_t N, Vec3 nx, Vec3 ny, Vec3 x0) {
  /*
    Returns a triangular mesh of a diamond with N^2 points.
    Each edge of the diamond has length 1.
    [out] m:  mesh of the diamond
    [in]  N:  number of points on an edge
    [in]  nx: unit normal, direction for x
    [in]  ny: unit normal, direction for y
    [in]  x0: bottom left starting point
  */
  
  // Build the vertices
  m.positions.resize(N*N);
  float step_size = 1.0/(N-1);
  for(int i=0; i<N; ++i)
    for(int j=0; j<N; ++j)
      m.positions[i+N*j] = x0 + step_size * (float(i)*nx + float(j)*ny);

  // Build the triangles
  m.indices.resize(3*2*(N-1)*(N-1));
  int try_id = 0;
  for(int i=0; i<N-1; ++i)
    for(int j=0; j<N-1; ++j) {
      int id_v_bl = i+N*j; // bottom left vtx
      
      m.indices[3*try_id] = id_v_bl;
      m.indices[3*try_id+1] = id_v_bl+1;
      m.indices[3*try_id+2] = id_v_bl+1+N;
      ++try_id;
      
      m.indices[3*try_id] = id_v_bl;
      m.indices[3*try_id+1] = id_v_bl+N;
      m.indices[3*try_id+2] = id_v_bl+N+1;
      ++try_id;
    }

  /* Assign types to the vertices */
  
  // Type 0 = inside diamond
  m.vertices_of_type[0].resize((N-2)*(N-2));
  int id_type_0 = 0;
  for(int i=1; i<N-1; ++i)
    for(int j=1; j<N-1; ++j) 
      m.vertices_of_type[0][id_type_0++] = i+N*j;

  //Type 1 = on the interior of an outer edge of the diamond
  m.vertices_of_type[1].resize(4*(N-2));
  int id_type_1 = 0;
  for(int i=1; i<N-1; ++i) {
    m.vertices_of_type[1][id_type_1++] = i;
    m.vertices_of_type[1][id_type_1++] = i+N*(N-1);
  }
  for(int j=1; j<N-1; ++j) {
    m.vertices_of_type[1][id_type_1++] = N*j;
    m.vertices_of_type[1][id_type_1++] = N-1+N*j;
  }

  // Type 2 = extremal vertex of the diamond
  m.vertices_of_type[2].resize(4);
  m.vertices_of_type[2][0] = 0;
  m.vertices_of_type[2][1] = N-1;
  m.vertices_of_type[2][2] = N*(N-1);
  m.vertices_of_type[2][3] = (N-1)+N*(N-1);
};

void load_cube(Mesh &m, size_t subdiv) {
  /*
    Returns a triangular mesh of the unit cube.
    [out] m:      mesh of the diamond
    [in]  subdiv: number of points on an edge of the cube
  */

  /* --- Construct each face individually --- */
  size_t N = subdiv;
  _face_mesh faces[6];
  Vec3 nx[6] = {
    Vec3::XAxis,
    Vec3::YAxis,
    Vec3::XAxis,
    -Vec3::XAxis,
    -Vec3::YAxis,
    -Vec3::XAxis,
  };
  Vec3 ny[6] = {
    Vec3::ZAxis,
    Vec3::ZAxis,
    Vec3::YAxis,
    -Vec3::ZAxis,
    -Vec3::ZAxis,
    -Vec3::YAxis,
  };
  Vec3 x0[6] = {
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {1, 1, 1},
    {1, 1, 1},
    {1, 1, 1},
  };
  for(int face_id=0; face_id<6; ++face_id)
    _load_diamond(faces[face_id], N, nx[face_id], ny[face_id], x0[face_id]);

  
  /* ----- Create the mesh of the cube ----- */
  
  /* --- Vertices update --- */
  size_t tot_nb_int_vtx = 6*(N-2)*(N-2);
  size_t tot_nb_edge_vtx = 12*(N-2);
  size_t tot_nb_extr_vtx = 8;
  size_t tot_nb_vtx = tot_nb_int_vtx + tot_nb_edge_vtx + tot_nb_extr_vtx;
  m.positions.resize(tot_nb_vtx);
  int cube_vtx_id = 0;
  float tol_sq = (0.5 / (N-1)) * (0.5 / (N-1));
  
  // (old_face_id, old_local_vtx_id) -> new_cube_vtx_id
  std::map<std::pair<int, int>, int> o2n_vtx;

  /* Update the vertices of type 0 (face_interior)
     --> no need to check for duplicates for them */
  // Number of vertices of type 0 per face
  int nb_vtx_type_0 = (N-2)*(N-2);
  for(int f=0; f<6; ++f)
    for(int v_t0_id=0; v_t0_id<nb_vtx_type_0; ++v_t0_id) {
      int v = faces[f].vertices_of_type[0][v_t0_id];
      for(int k=0; k<3; ++k)
	m.positions[cube_vtx_id][k] = faces[f].positions[v][k];
      o2n_vtx.insert({{f, v}, cube_vtx_id});
      ++cube_vtx_id;
    }

  /* --- Update the vertices of type 1 (outer edge interior) and 2 (face extremal vertices) ---*/
  // Number of vertices f type t-1 per face
  size_t nb_t_vtx[2] = {4*(N-2), 4};
  TArray<bool> vtx_was_handled[2][6];
  for(int i=0; i<2; ++i) {
    for(int f=0; f<6; ++f) {
      vtx_was_handled[i][f].resize(nb_t_vtx[i]);
      for(int k=0; k<nb_t_vtx[i]; ++k) {
	vtx_was_handled[i][f][k] = false;
      }
    }
  }
  
  for(int type=1; type<3; ++type) {
    // Number of duplicates for each vertex of current type
    int nb_eq_vtx = type+1;
      
    for(int f1=0; f1<6; ++f1)
      for(int t_id_1=0; t_id_1<nb_t_vtx[type-1]; ++t_id_1) {
	int v1 = faces[f1].vertices_of_type[type][t_id_1];

	// If v1 not handled, then add its coordinates to the cube mesh
	if (!vtx_was_handled[type-1][f1][t_id_1]) {
	  for(int k=0; k<3; ++k)
	    m.positions[cube_vtx_id][k] = faces[f1].positions[v1][k];
	  vtx_was_handled[type-1][f1][t_id_1] = true;
	  o2n_vtx.insert({{f1, v1}, cube_vtx_id});
	  ++cube_vtx_id;
	
	  // Find the duplicates of v1 and handle them
	  int nb_dup_found = 0;
	  bool all_dup_found = false;
	  for(int f2=f1+1; f2<6 && !all_dup_found; ++f2)
	    for(int t_id_2=0; t_id_2<nb_t_vtx[type-1] && !all_dup_found; ++t_id_2) {
	      int v2 = faces[f2].vertices_of_type[type][t_id_2];
	  
	      float dist_sq = norm2(faces[f1].positions[v1] - faces[f2].positions[v2]);
	      if (dist_sq < tol_sq) {
		vtx_was_handled[type-1][f2][t_id_2] = true;
		o2n_vtx.insert({{f2, v2}, cube_vtx_id-1});

		// If all duplicates found then end the search
		++nb_dup_found;
		all_dup_found = (nb_dup_found >= nb_eq_vtx-1);
	      }
	    }
	}
      }
  }
  
  assert(cube_vtx_id == tot_nb_vtx);
  

  /* --- Triangles update --- */
  size_t nb_tri_face = 2*(N-1)*(N-1);
  size_t tot_nb_tri = 6*2*(N-1)*(N-1);
  m.indices.resize(3*6*nb_tri_face);
  int cube_indices_id = 0;
  
  for(int f=0; f<6; ++f)
    for(int loc_tri=0; loc_tri<nb_tri_face; ++loc_tri)
      for(int k=0; k<3; ++k) {
	int loc_v = faces[f].indices[3*loc_tri+k];
	auto it = o2n_vtx.find({f, loc_v});
	assert(it != o2n_vtx.end());
	m.indices[cube_indices_id++] = it->second;
      }
};
