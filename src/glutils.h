#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <cstdlib>

#include <iostream>


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
GLuint uploadVertices(unsigned int size); // Uploads no data

// Compiles the sources and links them to a program, whose reference it returns
GLuint createShaderProgram(const GLchar* vertSrc, const GLchar* fragSrc, GLuint& vertexShaderRef, GLuint& fragmentShaderRef);

GLuint createVertexArrayObject();

void setVertexAttribPointer(GLuint shaderProgram, const char* name, GLint numComponents, GLenum type, GLsizei stride, int offset);


/* Does not bind VAO or buffer! */
void setFormat(const char *format, GLuint shaderProgram);

/* BUFFER WRITER */
// Helps write to a mapped buffer
// WARNING: Map the correct buffer to GL_ARRAY_BUFFER first
template <typename T>
class BufferWriter
{
public:
    // size: Size of buffer, in indices.
    BufferWriter(int size)
    {
        unmapped = false;
		ptr = static_cast<T*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, size * sizeof(T), GL_MAP_WRITE_BIT));
        if (ptr == nullptr) assert(!"Couldn't map buffer range - possible too large range.");
        this->size = 0;
        max_size = size;
    }
    BufferWriter(const BufferWriter& copy) = delete;

    ~BufferWriter() {
        assert(unmapped);
    }
    // Write to buffer
    void write(T a) {
        assert(!unmapped);
        if (size > max_size - 1) {
            assert(!"Buffer overflow.");
            exit(1);
            return;
        }
        *ptr = a; ptr ++;   
        size ++;
    }
    void write(T a, T b) {
        assert(!unmapped);
        if (size > max_size - 2) {
            assert(!"Buffer overflow.");
            return;
        }
        // std::cout << "WRITE " << size << " vs " << max_size << std::endl;
        *(ptr ++) = a; *(ptr ++) = b;
        size += 2;
    }
    void write(T a, T b, T c) {
        assert(!unmapped);
        if (size > max_size - 3) {
            assert(!"Buffer overflow.");
            exit(1);
            return;
        }
        *(ptr ++) = a; *(ptr ++) = b; *(ptr ++) = c;
        size += 3;
    }

    T* get_ptr() {
        assert(!unmapped);
        return ptr;
    }
    int get_current_size() {
        return size;
    }

    void unmap() {
        unmapped = true;
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }
private:
    bool unmapped;
    T* ptr;
    int size;
    int max_size;
};
