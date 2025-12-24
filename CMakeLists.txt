cmake_minimum_required(VERSION 3.25)

project(MU5MAM30 CXX)
set(CMAKE_BUILD_TYPE Debug)


add_executable(main)
target_sources(main
	PRIVATE src/main.cpp)
target_link_libraries(main
	PRIVATE
	mesh
	sphere
	array
	P1
)

add_library(mesh)
target_sources(mesh
        PRIVATE
        src/mesh/mesh.cpp

        PUBLIC
        FILE_SET HEADERS
	BASE_DIRS include/mesh
	FILES include/mesh/mesh.h
        )
target_link_libraries(mesh
	PRIVATE array)
add_library(sphere)
target_sources(sphere
        PRIVATE
        src/mesh/sphere.cpp

        PUBLIC
        FILE_SET HEADERS
	BASE_DIRS include/mesh
	FILES include/mesh/sphere.h
        )
target_link_libraries(sphere
	PRIVATE mesh array)
add_library(P1)
target_sources(P1
        PRIVATE
        src/fem/P1.cpp

        PUBLIC
        FILE_SET HEADERS
	BASE_DIRS include/fem
	FILES include/fem/P1.h
        )
target_link_libraries(P1
	PRIVATE fem_matrix sparse_matrix matrix array mesh)
add_library(fem_matrix)
target_sources(fem_matrix
        PRIVATE
        src/matrix/fem_matrix.cpp

        PUBLIC
        FILE_SET HEADERS
	BASE_DIRS include/matrix
	FILES include/matrix/fem_matrix.h
        )
target_link_libraries(fem_matrix
	PRIVATE array matrix mesh)
add_library(sparse_matrix)
target_sources(sparse_matrix
        PRIVATE
        src/matrix/sparse_matrix.cpp

        PUBLIC
        FILE_SET HEADERS
	BASE_DIRS include/matrix/
	FILES include/matrix/sparse_matrix.h
        )
target_link_libraries(sparse_matrix
	PRIVATE array stb_image_write)
add_library(array INTERFACE)
target_include_directories(array INTERFACE include/common)
target_sources(array
	INTERFACE
        FILE_SET HEADERS
	BASE_DIRS include/common
	FILES include/common/array.h
        )
add_library(matrix INTERFACE)
target_include_directories(array INTERFACE include/matrix)
target_sources(array
	INTERFACE
        FILE_SET HEADERS
	BASE_DIRS include/matrix
	FILES include/matrix/matrix.h
        )
add_library(stb_image_write INTERFACE)
target_include_directories(stb_image_write INTERFACE extern)
target_sources(stb_image_write
	INTERFACE
        FILE_SET HEADERS
	BASE_DIRS extern
	FILES extern/stb_image_write.h
        )
