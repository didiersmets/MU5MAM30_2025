#pragma once

/******************************************************************************
 *
 * OPENGL HEADERS & EXTENSION  HEADERS
 *
 *****************************************************************************/

#if defined __linux__
	#ifndef GL_GLEXT_PROTOTYPES
		#define GL_GLEXT_PROTOTYPES
	#endif
	#include <GL/gl.h>
	#include <GL/glext.h>
#elif defined __APPLE__
	#define GL_SILENCE_DEPRECATION 1
	#include <OpenGL/gl3.h>
	#include <OpenGL/glext.h>
#elif defined _WIN32
	#include "glad/glad.h" // utilisation de glad pour windows
#else
	#error "This platform doesn't seem supported yet."
#endif
