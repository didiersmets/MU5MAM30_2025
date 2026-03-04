#pragma once
#define USE_P2 true

#include <stdint.h>

#include "array.h"
#include "vec3.h"

#if USE_P2
#include <unordered_map>
#endif

struct Mesh {
  TArray<Vec3> positions;
  TArray<uint32_t> indices;
  TArray<float> attr;
#if USE_P2
  std::unordered_map<uint64_t, uint32_t> edge2dof;
  TArray<uint32_t> e2vtx;
#endif
  size_t vertex_count() const { return positions.size; }
  size_t index_count() const { return indices.size; }
  size_t triangle_count() const { return indices.size / 3; }
};

static inline uint64_t pack(uint32_t a, uint32_t b) noexcept {
  return (uint64_t(a) << 32) | uint64_t(b);
}
static inline void unpack(uint64_t p, uint32_t &a, uint32_t &b) noexcept {
  a = uint32_t(p >> 32);
  b = uint32_t(p & 0xFFFFFFFFu);
}
