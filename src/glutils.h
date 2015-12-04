#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>


void GLFW_boilerPlate(GLFWwindow **window, GLFWerrorfun error_callback);
/* Creates view matrix which really only offsets. */
glm::mat4 viewMatrix2D(float centerX, float centerY, float scaleX, float scaleY);

/* Assumes that the viewing area is symmetric (right = - left, top = - bottom) */
glm::mat4 ortho2D(float width, float height, float far, float near);

GLuint uploadElements(GLuint elements[], unsigned int size);

template <typename T>
GLuint uploadVertices(T data[], unsigned int size)
{
	GLuint reference;
	glGenBuffers(1,&reference);
	glBindBuffer(GL_ARRAY_BUFFER, reference);
	glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW);
	return reference;
}

// Compiles the sources and links them to a program, whose reference it returns
GLuint createShaderProgram(const GLchar* vertSrc, const GLchar* fragSrc, GLuint& vertexShaderRef, GLuint& fragmentShaderRef);

GLuint createVertexArrayObject();

void setVertexAttribPointer(GLuint shaderProgram, const char* name, GLint numComponents, GLenum type, GLsizei stride, int offset);


/* Does not bind VAO or buffer! */
void setFormat(const char *format, GLuint shaderProgram);

// Helps write to a mapped buffer:
// TODO Doesn't unmap...
template <typename T>
class BufferWriter
{
public:
    // size: Size of buffer, in indices.
    BufferWriter(int size)
    {
		ptr = static_cast<T*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, size * sizeof(T), GL_MAP_WRITE_BIT));
    }
    // Write to buffer
    void write(T a) {
        *ptr = a; ptr ++;   
    }
    void write(T a, T b) {
        *(ptr ++) = a; *(ptr ++) = b;
    }
    void write(T a, T b, T c) {
        *(ptr ++) = a; *(ptr ++) = b; *(ptr ++) = c;
    }
private:
    T* ptr;
};
