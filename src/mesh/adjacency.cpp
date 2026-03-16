#include "mesh.h"
#include "adjacency.h"

// VTAdjacency::VTAdjacency(const Mesh &m)
// {
// 	/* Your implementation goes here */
// }


VTAdjacency::VTAdjacency(const Mesh &m)
{
	size_t vtx_count = m.vertex_count(); //顶点数量
	size_t tri_count = m.triangle_count(); //三角形数量
	size_t idx_count = 3 * tri_count; //三角形所有顶点索引数量
	degree.resize(vtx_count); //resize来自array.h文件,用于调整数组的大小为指定值 
	offset.resize(vtx_count);
	vtri.resize(idx_count);

	for (size_t i = 0; i < vtx_count; ++i) {
		degree[i] = 0;
		offset[i] = 0;
	}
	// 遍历所有顶点索引,每当一个顶点在索引出现,degree就+1,degree[i]表示顶点i所连接的三角形数量
	// 俩三角形的m.indices=[0,1,2,1,2,3] -> degree[0]=2 顶点0涉及俩三角形 
	for (size_t i = 0; i < idx_count; ++i) {
		degree[m.indices[i]]++;
	}

	//累加前一个顶点的degree,计算每个顶点在vtri中的起始位置
	for (size_t v = 1; v < vtx_count; ++v) {
		offset[v] = offset[v - 1] + degree[v - 1];
	}

	// 构建vtri表(顶点和相邻点信息)
	// 遍历每个三角形,根据其三个顶点,vtri存储相邻的两个顶点
	for (size_t t = 0; t < tri_count; ++t) {
		uint32_t a = m.indices[3 * t + 0];
		uint32_t b = m.indices[3 * t + 1];
		uint32_t c = m.indices[3 * t + 2];
		vtri[offset[a]++] = { b, c };
		vtri[offset[b]++] = { c, a };
		vtri[offset[c]++] = { a, b };
	}

	//还原offset为每个顶点的起始位置
	for (size_t v = 0; v < vtx_count; ++v) {
		offset[v] -= degree[v];
	}
	assert(offset[vtx_count - 1] + degree[vtx_count - 1] == idx_count);//验证
}

/* vtri = [
    { next: 1, prev: 2 }, // 顶点 0
    { next: 2, prev: 0 }, // 顶点 1 (第一次出现)
    { next: 3, prev: 2 }, // 顶点 1 (第二次出现)
    { next: 0, prev: 1 }, // 顶点 2 (第一次出现)
    { next: 1, prev: 3 }, // 顶点 2 (第二次出现)
    { next: 2, prev: 1 }, // 顶点 3
	];
 */