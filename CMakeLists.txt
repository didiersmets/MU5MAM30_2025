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
	# NOTE : If you want to test your code step by step, before completion
	#        of all the files with "holes" (i.e. sections labeled 
	#        /* Your implementation goes here */ ), you can simply prefix
	#        the corresponding unfinished source files below with the
	#        solutions directory. 
	#        If P1.cpp has not been completed yet e.g., you simply replace :
	#	    src/fem/P1.cpp
	#        by
	#           solutions/src/fem/P1.cpp
	#        The files mass.h and stiffness.h are specific since they are
	#        header files, you should complete (or copy them from solutions)
	#        before building the project.
	src/fem/P1.cpp
	src/fem/P2.cpp #ajoute de P2 Aurélien
	src/fem/poisson.cpp
	src/fem/navier_stokes.cpp
	src/linalg/conjugate_gradient.cpp
	src/linalg/cholesky.cpp #Ajout de cholesky
	src/matrix/fem_matrix.cpp
	src/matrix/sparse_matrix.cpp
	#src/matrix/spy.cpp
	src/mesh/adjacency.cpp
	src/mesh/cube.cpp
	src/mesh/cube2.cpp
	src/mesh/sphere.cpp
	src/mesh/duplicate_verts.cpp
	src/mesh/mesh_io.cpp
	src/mesh/mesh_bounds.cpp
    src/mesh/adjacency_P2.cpp
	include/mesh/adjacency_P2.h
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
	-std=c++17 -Wall -Wextra -march=native -Wno-unknown-pragmas
	-Wno-missing-field-initializers
)
if (USE_OPENMP AND OpenMP_CXX_FOUND)
	target_compile_options(nanofem PRIVATE -DUSE_OPENMP)
	target_link_libraries(nanofem PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(nanofem PUBLIC m)
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
	src/viewer/glad.c #rajouté pour windows par Aurélien

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
add_executable(test_NS
	src/common/logging.cpp
	src/bin/test_navier_stokes.cpp
)
target_compile_options(test_NS PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_NS PRIVATE nanofem nanoviewer tiny_expr)

###############################################################################
add_executable(test_cholesky
	src/bin/test_cholesky.cpp
)
target_compile_options(test_cholesky PRIVATE 
	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
)
target_link_libraries(test_cholesky PRIVATE nanofem)
###############################################################################

#
# copy the shaders directory to the build tree so that the hard-coded
# relative paths used in the viewer programs ("./shaders/fem.vert" etc.)
# actually resolve when the executables are launched from the build
# directory.  Without this CMake command the project root is the only
# place that contains the "shaders" folder which explains the runtime
# error shown in the bug report:
#
#   ERROR: Vertex shader ./shaders/fem.vert not found!
#
# Two alternatives are possible: either change the sources to compute a
# full path to the data from argv[0] / executable location, or make the
# shaders available in the working directory.  the latter is simpler for
# a student project so we copy the directory during the build.
#
set(SHADER_SRC_DIR ${CMAKE_SOURCE_DIR}/shaders)

add_custom_target(copy_shaders ALL
	COMMAND ${CMAKE_COMMAND} -E copy_directory
			"${SHADER_SRC_DIR}" "${CMAKE_BINARY_DIR}/shaders"
	COMMENT "copying shader assets to build directory"
)

add_dependencies(test_NS copy_shaders)
add_dependencies(test_poisson copy_shaders)
###############################################################################

###############################################################################
#add_executable(test_ordering
#	src/common/logging.cpp
#	src/common/chrono.cpp
#	src/bin/test_ordering.cpp
#)
#target_compile_options(test_ordering PRIVATE 
#	-std=c++17 -Wall -Wextra -Wno-missing-field-initializers
#)
#target_link_libraries(test_ordering PRIVATE nanofem)
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

