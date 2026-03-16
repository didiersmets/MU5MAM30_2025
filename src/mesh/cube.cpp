#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "cube.h"
#include "duplicate_verts.h"
#include "math_utils.h"
#include "mesh.h"
#include "sys_utils.h"
#include "vec3.h"

static void load_cube_vertices(Vec3 *pos, size_t subdiv);
static void load_cube_indices(uint32_t *idx, size_t subdiv);

int load_cube(Mesh &m, size_t subdiv)
{
	/* Check subdiv is reasonable and return error if not */
	if (subdiv <= 0 || subdiv > (1 << 14) /* 16K */) {
		return (-1);
	}

	size_t n = subdiv + 1;

	/* Reserve memory for vertices and indices */
	m.positions.resize(6 * POW2(n));
	m.indices.resize(36 * POW2(subdiv));

	/* First build vertices as six unattached faces of n^2 vertices each */
	/* See below for implementation */
	load_cube_vertices(m.positions.data, subdiv);

	/* Build corresponding triangulation indices */
	/* See below for implementation */
	load_cube_indices(m.indices.data, subdiv);

	/* Finally attach faces between themselves */
	/* Implementation in src/duplicate_verts.cpp */
	remove_duplicate_vertices(m);

	return (0);
}

// static void load_cube_vertices(Vec3 *pos, size_t subdiv)
// {
// 	/* Your implementation goes here */
// }

// static void load_cube_indices(uint32_t *idx, size_t subdiv)
// {
// 	/* Your implementation goes here */
// }


// 生成立方体的顶点
// pos存储生成的顶点的数组	subdiv决定每个面的顶点数量
static void load_cube_vertices(Vec3 *pos, size_t subdiv)
{
	size_t n = subdiv + 1; //每个面在每个方向上的顶点数量

	/* First build vertices as six unattached faces of n^2 vertices each */
	size_t voff[6]; //每个面的顶点在pos中的起始偏移量(共6个面)
	for (int f = 0; f < 6; ++f) {
		voff[f] = f * POW2(n);
	}

	float *dir = (float *)safe_malloc(n * sizeof(float)); //顶点坐标沿轴的正向分布
	float *rev = (float *)safe_malloc(n * sizeof(float)); //顶点坐标沿轴的反向分布 (将dir反过来就是rev)
	for (size_t i = 0; i < n; ++i) {
		dir[i] = (2 * (float)i - subdiv) / subdiv;
	}
	for (size_t i = 0; i < n; ++i) {
		rev[i] = dir[subdiv - i];
	}

	for (size_t y = 0; y < n; ++y) {
		for (size_t x = 0; x < n; ++x) {
			pos[voff[0]++] = { dir[x], -1, dir[y] }; /* Front 存在pos[0]-pos[15]中  */
			pos[voff[1]++] = { rev[x], +1, dir[y] }; /* Back   */
			pos[voff[2]++] = { -1, rev[x], dir[y] }; /* Left   */
			pos[voff[3]++] = { +1, dir[x], dir[y] }; /* Right  */
			pos[voff[4]++] = { dir[x], rev[y], -1 }; /* Bottom */
			pos[voff[5]++] = { dir[x], dir[y], +1 }; /* Top    */
		}
	}
	free(dir);
	free(rev);
}

//生成立方体的三角形索引
static void load_cube_indices(uint32_t *idx, size_t subdiv)
{
	size_t n = subdiv + 1;
	/* Build corresponding triangulation indices */
	for (int f = 0; f < 6; f++) { //遍历6个面
		size_t offset = f * POW2(n); //每个面在顶点数组中的起始偏移量
		for (size_t i = 0; i < subdiv; ++i) { 
			for (size_t j = 0; j < subdiv; ++j) {
				uint32_t base = (uint32_t)(i * n + j + offset); //当前网格左上角
				/* First tri 右上三角形 */
				*idx++ = base; //左上顶点
				*idx++ = base + 1;//右上顶点
				*idx++ = base + 1 + n; //右下顶点
				/* Second tri 左下三角形*/
				*idx++ = base;
				*idx++ = base + 1 + n;
				*idx++ = base + n; //左下顶点
			}
		}
	}
}