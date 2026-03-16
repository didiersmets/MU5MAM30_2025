cmake_minimum_required(VERSION 3.14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

project(MU5MAM30 LANGUAGES CXX C)

find_package(OpenGL REQUIRED)

###############################################################################
set(USE_OPENMP FALSE)
if (USE_OPENMP)
	find_package(OpenMP)
endif()
add_library(nanofem
	src/fem/P1.cpp
	src/fem/poisson.cpp
	src/linalg/conjugate_gradient.cpp
	src/linalg/projected_conjugate_gradient.cpp
	src/matrix/fem_matrix.cpp
	src/matrix/sparse_matrix.cpp
	src/mesh/adjacency.cpp
	src/mesh/cube.cpp
	src/mesh/cube2.cpp
	src/mesh/sphere.cpp
	src/mesh/half_sphere.cpp
	src/mesh/duplicate_verts.cpp
	src/mesh/mesh_io.cpp
	src/mesh/mesh_bounds.cpp
)
target_include_directories(nanofem PUBLIC 
	include/common 
	include/fem
	include/linalg
	include/matrix 
	include/mesh 
	extern
)
target_compile_options(nanofem PRIVATE 
	-std=c++17 -Wall -Wextra -ffast-math -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
if (USE_OPENMP AND OpenMP_CXX_FOUND)
	target_compile_options(nanofem PRIVATE -DUSE_OPENMP)
	target_link_libraries(nanofem PUBLIC OpenMP::OpenMP_CXX)
endif()

# CG variant (uses conjugate_gradient + navier_stokes_CG)
add_library(nanofem_CG
	src/fem/navier_stokes_CG.cpp
	src/linalg/conjugate_gradient.cpp
)
target_include_directories(nanofem_CG PUBLIC 
	include/common 
	include/fem
	include/linalg
	include/matrix 
	include/mesh 
	extern
)
target_compile_options(nanofem_CG PRIVATE 
	-std=c++17 -Wall -Wextra -ffast-math -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
target_link_libraries(nanofem_CG PUBLIC nanofem)

# CG non-penetrating BC variant
add_library(nanofem_CG_nonpenetrating
	src/fem/navier_stokes_CG_nonpenetrating.cpp
	src/linalg/conjugate_gradient.cpp
	src/linalg/projected_conjugate_gradient.cpp
)
target_include_directories(nanofem_CG_nonpenetrating PUBLIC 
	include/common 
	include/fem
	include/linalg
	include/matrix 
	include/mesh 
	extern
)
target_compile_options(nanofem_CG_nonpenetrating PRIVATE 
	-std=c++17 -Wall -Wextra -ffast-math -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
target_link_libraries(nanofem_CG_nonpenetrating PUBLIC nanofem)

# CG no-slip BC variant (non-penetrating + no-slip)
add_library(nanofem_CG_noslip
	src/fem/navier_stokes_CG_noslip.cpp
	src/linalg/conjugate_gradient.cpp
	src/linalg/projected_conjugate_gradient.cpp
)
target_include_directories(nanofem_CG_noslip PUBLIC 
	include/common 
	include/fem
	include/linalg
	include/matrix 
	include/mesh 
	extern
)
target_compile_options(nanofem_CG_noslip PRIVATE 
	-std=c++17 -Wall -Wextra -ffast-math -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
target_link_libraries(nanofem_CG_noslip PUBLIC nanofem)

# Cholesky variant (uses cholesky + navier_stokes_cholesky)
add_library(nanofem_cholesky
	src/fem/navier_stokes_cholesky.cpp
	src/linalg/cholesky.cpp
	src/linalg/conjugate_gradient.cpp
)
target_include_directories(nanofem_cholesky PUBLIC 
	include/common 
	include/fem
	include/linalg
	include/matrix 
	include/mesh 
	extern
)
target_compile_options(nanofem_cholesky PRIVATE 
	-std=c++17 -Wall -Wextra -ffast-math -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
target_link_libraries(nanofem_cholesky PUBLIC nanofem)
###############################################################################

###############################################################################
add_library(nanoviewer
	src/mesh/mesh_gpu.cpp
	src/viewer/viewer.cpp
	src/viewer/camera.cpp
	src/viewer/frustum.cpp
	src/viewer/mouse.cpp
	src/viewer/trackball.cpp
	src/viewer/shaders.cpp
)
target_include_directories(nanoviewer PUBLIC 
	include/common 
	include/mesh
	include/viewer
	extern
)
target_compile_options(nanoviewer PRIVATE
	-std=c++17 -Wall -Wextra
)
target_link_libraries(nanoviewer PUBLIC
	glfw
	imgui_glfw_opengl3
	OpenGL::GL
)
###############################################################################

###############################################################################
add_executable(test_poisson
	src/common/logging.cpp
	src/bin/test_poisson.cpp
)
target_compile_options(test_poisson PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_poisson PRIVATE nanofem nanoviewer tiny_expr) 
###############################################################################

###############################################################################
add_library(nanofem_poisson_dirichlet
	src/fem/poisson_non_homogene_Dirichlet.cpp
	src/linalg/projected_conjugate_gradient.cpp
)
target_include_directories(nanofem_poisson_dirichlet PUBLIC
	include/common
	include/fem
	include/linalg
	include/matrix
	include/mesh
	extern
)
target_compile_options(nanofem_poisson_dirichlet PRIVATE
	-std=c++17 -Wall -Wextra -ffast-math -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
target_link_libraries(nanofem_poisson_dirichlet PUBLIC nanofem)
###############################################################################

###############################################################################
add_executable(test_poisson_dirichlet
	src/common/logging.cpp
	src/bin/test_poisson_dirichlet.cpp
)
target_compile_options(test_poisson_dirichlet PRIVATE
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_poisson_dirichlet PRIVATE nanofem nanofem_poisson_dirichlet nanoviewer tiny_expr)
###############################################################################

###############################################################################
add_executable(test_NS_CG
	src/common/logging.cpp
	src/bin/test_navier_stokes_CG.cpp
)
target_compile_options(test_NS_CG PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_NS_CG PRIVATE nanofem nanofem_CG nanoviewer tiny_expr)
###############################################################################

###############################################################################
add_executable(test_NS_CG_nonpenetrating
	src/common/logging.cpp
	src/bin/test_navier_stokes_CG_nonpenetrating.cpp
)
target_compile_options(test_NS_CG_nonpenetrating PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_NS_CG_nonpenetrating PRIVATE nanofem nanofem_CG_nonpenetrating nanoviewer tiny_expr)
###############################################################################

###############################################################################
add_executable(test_NS_CG_noslip
	src/common/logging.cpp
	src/bin/test_navier_stokes_CG_noslip.cpp
)
target_compile_options(test_NS_CG_noslip PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_NS_CG_noslip PRIVATE nanofem nanofem_CG_noslip nanoviewer tiny_expr)
###############################################################################

###############################################################################
add_executable(test_NS_cholesky
	src/common/logging.cpp
	src/bin/test_navier_stokes_cholesky.cpp
)
target_compile_options(test_NS_cholesky PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_NS_cholesky PRIVATE nanofem nanofem_cholesky nanoviewer tiny_expr)
###############################################################################

###############################################################################
add_library(imgui_glfw_opengl3
	extern/imgui/imgui.cpp
	extern/imgui/imgui_draw.cpp
	extern/imgui/imgui_tables.cpp
	extern/imgui/imgui_widgets.cpp
	extern/imgui/imgui_impl_glfw.cpp
	extern/imgui/imgui_impl_opengl3.cpp
)
target_include_directories(imgui_glfw_opengl3 PUBLIC
	extern/glfw/include
)
###############################################################################

###############################################################################
add_library(tiny_expr
	extern/tiny_expr/tinyexpr.c
)
target_include_directories(tiny_expr PUBLIC
	extern
)
###############################################################################

add_subdirectory(extern/glfw)