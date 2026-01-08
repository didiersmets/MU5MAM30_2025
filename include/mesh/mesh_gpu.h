#pragma once

#include "common/gl_utils.h"
#include "mesh/mesh.h"

struct GPUMesh
{
  const Mesh* m;
  GLuint      pos_vbo;
  GLuint      idx_vbo;
  GLuint      attr_vbo;
  GLuint      vao;
  void        upload();
  void        update_attr();
  void        draw() const;
};
