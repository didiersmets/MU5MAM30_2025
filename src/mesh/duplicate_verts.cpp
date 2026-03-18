#include "duplicate_verts.h"
#include <algorithm> // For std::sort
#include <vector>    // For std::vector (used internally for sorting)

#include "mesh.h"
#include "vec3.h"

// Helper struct to store an index and verify which position it points to
struct VertRef {
    const Vec3* posArray; // Pointer to the start of the positions array
    uint32_t index;       // The index of this specific vertex

    // Operator to compare two vertices lexicographically (X, then Y, then Z)
    bool operator<(const VertRef& other) const {
        const Vec3& a = posArray[index];
        const Vec3& b = other.posArray[other.index];
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    }
    
    // Check for equality
    bool operator==(const VertRef& other) const {
        const Vec3& a = posArray[index];
        const Vec3& b = other.posArray[other.index];
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

size_t build_position_remap(const Vec3 *pos, size_t count, uint32_t *remap) {
    // 1. Create a list of references to all vertices
    std::vector<VertRef> refs(count);
    for (uint32_t i = 0; i < count; ++i) {
        refs[i].posArray = pos;
        refs[i].index = i;
    }

    // 2. Sort the references based on the XYZ position
    // This groups identical vertices together in the list
    std::sort(refs.begin(), refs.end());

    // 3. Build the remapping table
    // Iterate through the sorted list. If a vertex is the same as the previous one,
    // map it to the same original index.
    size_t unique_count = 0;
    for (size_t i = 0; i < count; ++i) {
        uint32_t original_idx = refs[i].index;

        if (i > 0 && refs[i] == refs[i-1]) {
            // This is a duplicate of the previous vertex in the sorted list.
            // Point it to the same remapped index as the previous one.
            uint32_t prev_original_idx = refs[i-1].index;
            remap[original_idx] = remap[prev_original_idx];
        } else {
            // This is a new, unique vertex (or the very first one).
            remap[original_idx] = original_idx;
            unique_count++;
        }
    }

    return unique_count;
}

void remove_duplicate_vertices(Mesh &m) {
    size_t old_count = m.vertex_count();
    if (old_count == 0) return;

    // 1. Calculate the remapping array
    std::vector<uint32_t> remap(old_count);
    
    build_position_remap(m.positions.data, old_count, remap.data());

    // 2. Create the new, compacted position array
    TArray<Vec3> new_positions;
    std::vector<uint32_t> old_to_new_dense(old_count); 
    
    for (size_t i = 0; i < old_count; ++i) {
        if (remap[i] == i) {
            //add unique vertext to the list
            old_to_new_dense[i] = (uint32_t)new_positions.size;
            new_positions.push_back(m.positions[i]);
        }
    }

    // 3. Update the indices to point to the new dense locations
    for (size_t i = 0; i < m.indices.size; ++i) {
        uint32_t old_idx = m.indices[i];
        uint32_t canonical_idx = remap[old_idx]; // Find the unique version of this vertex
        m.indices[i] = old_to_new_dense[canonical_idx]; // Find where that unique version lives now
    }

    // 4. Swap in the new positions
    auto old_data = m.positions.data; 
    
    m.positions.data = new_positions.data;
    m.positions.size = new_positions.size;

    // 3. Nullify the local array so it doesn't delete the data when it goes out of scope
    new_positions.data = nullptr;
    new_positions.size = 0;
   
}