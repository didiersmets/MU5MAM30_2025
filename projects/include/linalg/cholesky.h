#pragma once

#include <vector>
#include "sparse_matrix.h"


void elimination_tree(const CSRMatrix& A, std::vector<int>& parent);