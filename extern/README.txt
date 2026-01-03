This directory contains third party code that is used in some non core part of
the project : 

1) GLFW is an Open Source, multi-platform library for OpenGL (Graphics Library),
it is only used to make "user firendly" testing applications to vizualize meshes
and FEM solutions : it has functions to create windows, and interact with the 
mouse and the keyboard.

2) ImGUI stands for Immediate mode Graphical User Interface. It is used to
create controls in the graphical apps : buttons, sliders, checkboxes etc.

3) Stbi is a public domain library to easily load and write images in a number
of common formats (including png and jpeg). We shall only use it to vizualize
sparsity structures of FEM matrices depending on DOF ordering.

4) tiny_expr is single file math expression parser, it transforms a (valid)
math expression string into a computation tree at run time (i.e. live, with
no need of a compilation). It is used in graphical 
application tests in order to let the user choose its forcing term (i.e.
r.h.s. function f) or initial data (for evolution equations) in the form 
of a mathematical expression.

5) fast_obj is a single file loader for meshes in Wavefront OBJ 
(https://en.wikipedia.org/wiki/Wavefront_.obj_file), for these that
would like to test FEM codes on more fancy meshes than cubes of spheres.
Wavefront OBJ has been a popular mesh format for at least two decades.
 
