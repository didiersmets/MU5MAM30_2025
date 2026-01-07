#main:
#	g++ src/main.cpp src/mesh/mesh.cpp src/mesh/sphere.cpp -I./include -o build/main.exe

cm:
	rm -r build
	cmake -B build -DCMAKE_BUILD_TYPE=Debug

b :
	cmake --build build 
