#include "mesh.h"
#include <map>

namespace {
  struct face_mesh {
    TArray<Vec3> positions;
    TArray<uint32_t> indices;
    TArray<int> vertices_of_type[3];
  };
}


static void load_face_vertices(face_mesh &m,
			       size_t N,
			       const Vec3 &nx,
			       const Vec3 &ny,
			       const Vec3 &x0)
{
  m.positions.resize(N*N);
  float step_size = 1.0/(N-1);

  for(size_t i=0; i<N; ++i)
    for(size_t j=0; j<N; ++j)
      m.positions[i+N*j] = x0 + step_size * (float(i)*nx + float(j)*ny);
}

static void load_face_indices(face_mesh &m, size_t N) {
  m.indices.resize(3*2*(N-1)*(N-1));
  int try_id = 0;

  for(size_t i=0; i<N-1; ++i)
    for(size_t j=0; j<N-1; ++j) {
      /* Bottom left vertex index */
      int id_v_bl = i+N*j;

      m.indices[3*try_id] = id_v_bl;
      m.indices[3*try_id+1] = id_v_bl+1;
      m.indices[3*try_id+2] = id_v_bl+1+N;
      ++try_id;

      m.indices[3*try_id] = id_v_bl;
      m.indices[3*try_id+1] = id_v_bl+N;
      m.indices[3*try_id+2] = id_v_bl+N+1;
      ++try_id;
    }
}

static void load_face_vertices_types(face_mesh &m, size_t N) {
  /* Type 0 = inside face */
  m.vertices_of_type[0].resize((N-2)*(N-2));
  int id_type_0 = 0;
  for(size_t i=1; i<N-1; ++i)
    for(size_t j=1; j<N-1; ++j)
      m.vertices_of_type[0][id_type_0++] = i+N*j;

  /* Type 1 = on the interior of an outer edge of the face */
  m.vertices_of_type[1].resize(4*(N-2));
  int id_type_1 = 0;
  for(size_t i=1; i<N-1; ++i) {
    m.vertices_of_type[1][id_type_1++] = i;
    m.vertices_of_type[1][id_type_1++] = i+N*(N-1);
  }
  for(size_t j=1; j<N-1; ++j) {
    m.vertices_of_type[1][id_type_1++] = N*j;
    m.vertices_of_type[1][id_type_1++] = N-1+N*j;
  }

  /* Type 2 = extremal vertex of the face */
  m.vertices_of_type[2].resize(4);
  m.vertices_of_type[2][0] = 0;
  m.vertices_of_type[2][1] = N-1;
  m.vertices_of_type[2][2] = N*(N-1);
  m.vertices_of_type[2][3] = (N-1)+N*(N-1);
}

static void load_face(face_mesh &m,
		      size_t N,
		      const Vec3 &nx,
		      const Vec3 &ny,
		      const Vec3 &x0)
{
  /*
    Returns a triangular mesh of a face with N^2 points.
    Each edge of the face has length 1.
    [out] m:  mesh of the face
    [in]  N:  number of points on an edge
    [in]  nx: unit normal, direction for x
    [in]  ny: unit normal, direction for y
    [in]  x0: bottom left starting point
  */

  /* Build the vertices of the mesh of the face */
  load_face_vertices(m, N, nx, ny, x0);

  /* Build the triangles of the mesh of the face */
  load_face_indices(m, N);

  /* Assign a type to each vertex */
  load_face_vertices_types(m, N);
}

static void load_cube_vertices(TArray<Vec3> &pos,
			       size_t N,
			       const face_mesh faces[6],
			       std::map<std::pair<int, int>, int> &o2n_vtx)
{
  /* (done) Your implementation goes here */

  size_t tot_nb_int_vtx = 6*(N-2)*(N-2);
  size_t tot_nb_edge_vtx = 12*(N-2);
  size_t tot_nb_extr_vtx = 8;
  size_t tot_nb_vtx = tot_nb_int_vtx + tot_nb_edge_vtx + tot_nb_extr_vtx;
  pos.resize(tot_nb_vtx);
  int cube_vtx_id = 0;
  float tol_sq = (0.5 / (N-1)) * (0.5 / (N-1));

  /* Update the vertices of type 0 (face_interior)
     --> no need to check for duplicates for them */
  /* Number of vertices of type 0 per face */
  int nb_vtx_type_0 = (N-2)*(N-2);
  for(int f=0; f<6; ++f)
    for(int v_t0_id=0; v_t0_id<nb_vtx_type_0; ++v_t0_id) {
      int v = faces[f].vertices_of_type[0][v_t0_id];
      for(int k=0; k<3; ++k)
	pos[cube_vtx_id][k] = faces[f].positions[v][k];
      o2n_vtx.insert({{f, v}, cube_vtx_id});
      ++cube_vtx_id;
    }

  /* Update the vertices of type 1 (outer edge interior) and 2 (face extremal
     vertices) */
  /* Number of vertices f type t-1 per face */
  size_t nb_t_vtx[2] = {4*(N-2), 4};
  TArray<bool> vtx_was_handled[2][6];
  for(int i=0; i<2; ++i) {
    for(int f=0; f<6; ++f) {
      vtx_was_handled[i][f].resize(nb_t_vtx[i]);
      for(size_t k=0; k<nb_t_vtx[i]; ++k) {
	vtx_was_handled[i][f][k] = false;
      }
    }
  }

  for(int type=1; type<3; ++type) {
    /* Number of duplicates for each vertex of current type */
    int nb_eq_vtx = type+1;

    for(int f1=0; f1<6; ++f1)
      for(size_t t_id_1=0; t_id_1<nb_t_vtx[type-1]; ++t_id_1) {
	int v1 = faces[f1].vertices_of_type[type][t_id_1];

	/* If v1 not handled, then add its coordinates to the cube mesh */
	if (!vtx_was_handled[type-1][f1][t_id_1]) {
	  for(int k=0; k<3; ++k)
	    pos[cube_vtx_id][k] = faces[f1].positions[v1][k];
	  vtx_was_handled[type-1][f1][t_id_1] = true;
	  o2n_vtx.insert({{f1, v1}, cube_vtx_id});
	  ++cube_vtx_id;

	  /* Find the duplicates of v1 and handle them */
	  int nb_dup_found = 0;
	  bool all_dup_found = false;
	  for(int f2=f1+1; f2<6 && !all_dup_found; ++f2)
	    for(size_t t_id_2=0;
		t_id_2<nb_t_vtx[type-1] && !all_dup_found;
		++t_id_2) {
	      int v2 = faces[f2].vertices_of_type[type][t_id_2];

	      float dist_sq
		= norm2(faces[f1].positions[v1] - faces[f2].positions[v2]);
	      if (dist_sq < tol_sq) {
		vtx_was_handled[type-1][f2][t_id_2] = true;
		o2n_vtx.insert({{f2, v2}, cube_vtx_id-1});

		/* If all duplicates found then end the search */
		++nb_dup_found;
		all_dup_found = (nb_dup_found >= nb_eq_vtx-1);
	      }
	    }
	}
      }
  }

  assert(cube_vtx_id == (int)tot_nb_vtx);
}
static void load_cube_indices(TArray<uint32_t> &idx,
			      size_t N,
			      face_mesh faces[6],
			      std::map<std::pair<int, int>, int> &o2n_vtx)
{
  /* (done) Your implementation goes here */

  size_t nb_tri_face = 2*(N-1)*(N-1);
  idx.resize(3*6*nb_tri_face);
  int cube_indices_id = 0;

  for(int f=0; f<6; ++f)
    for(size_t loc_tri=0; loc_tri<nb_tri_face; ++loc_tri)
      for(int k=0; k<3; ++k) {
	int loc_v = faces[f].indices[3*loc_tri+k];
	auto it = o2n_vtx.find({f, loc_v});
	assert(it != o2n_vtx.end());
	idx[cube_indices_id++] = it->second;
      }
}


int load_cube(Mesh &m, size_t subdiv) {
  /*
    Returns a triangular mesh of the unit cube.
    [out] m:      mesh of the face
    [in]  subdiv: number of subdivision in each direction
  */

  /* Check subdiv is reasonable and return error if not */
  if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */) {
    return (-1);
  }

  size_t N = subdiv + 1;

  /* Construct each face individually */
  face_mesh faces[6];
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
    load_face(faces[face_id], N, nx[face_id], ny[face_id], x0[face_id]);


  /* --- Create the mesh of the cube --- */
  /* (old_face_id, old_local_vtx_id) -> new_cube_vtx_id connectivity */
  std::map<std::pair<int, int>, int> o2n_vtx;

  /* Update the vertices */
  load_cube_vertices(m.positions, N, faces, o2n_vtx);

  /* Update the indices */
  load_cube_indices(m.indices, N, faces, o2n_vtx);

  return (0);
};
