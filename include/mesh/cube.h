#pragma once
#include "mesh.h"

void rotate_y_to_x(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos);
void rotate_y_to_z(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos);
void flip_my_to_py(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos);
void flip_mx_to_px(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos);
void flip_mz_to_pz(const TArray<Vec3> &ref_pos, TArray<Vec3> &pos);

void load_cube(Mesh &m, size_t subdiv);
int load_cube_nested_dissect(Mesh &m, size_t subdiv);
