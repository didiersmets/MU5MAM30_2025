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
add_library(array INTERFACE)
target_include_directories(array INTERFACE include/common)
target_sources(array
	INTERFACE
        FILE_SET HEADERS
	BASE_DIRS include/common
	FILES include/common/array.h
        )
