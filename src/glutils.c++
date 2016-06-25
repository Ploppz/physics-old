#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <cassert>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <vector>

glm::mat4 viewMatrix2D(float centerX, float centerY, float scaleX, float scaleY)
{
	// data views the transpose of the actual matrix
	// glm interprets it in weird ways
	float matrixData[16] = {
		scaleX,	0, 0,	0,
		0,	scaleY, 0,	0,
		0,	0,	1,	0,
		-centerX,	-centerY,	0,	1
	};
	glm::mat4 matrix = glm::make_mat4x4(matrixData);
	return matrix;
}
glm::mat4 ortho2D(float width, float height, float far, float near)
{
	// Assumes that the viewing area is symmetric (right = - left, top = - bottom)
	// matrixData views the transpose of the actual matrix
	// glm interprets it in weird ways
	float matrixData[16] = {
		2/width, 0, 			0, 							0,
		0, 			 2/height,  0, 							0,
		0, 			 0,  			-2/(far - near), 			0,
		0, 			 0, 			-(far + near)/(far - near), 1
	};
	glm::mat4 matrix = glm::make_mat4x4(matrixData);
	return matrix;
}



typedef unsigned int uint;

GLuint uploadElements(GLuint elements[], unsigned int size)
{
	GLuint reference;
	glGenBuffers(1, &reference);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, reference);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, elements, GL_STATIC_DRAW);
	return reference;
}


GLuint uploadVertices(unsigned int size)
{
	GLuint reference;
	glGenBuffers(1,&reference);
	glBindBuffer(GL_ARRAY_BUFFER, reference);
	glBufferData(GL_ARRAY_BUFFER, size, NULL, GL_STREAM_DRAW);
	return reference;
}


// Compiles the sources and links them to a program, whose reference it returns

GLuint createShaderProgram(const GLchar* vertSrc, const GLchar* fragSrc, GLuint& vertexShaderRef, GLuint& fragmentShaderRef)
{	
	// COMPILING VERTEX SHADER
	vertexShaderRef = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShaderRef, 1, &vertSrc, NULL);
	glCompileShader(vertexShaderRef);

	// COMPILING FRAGMENT SHADER
	fragmentShaderRef = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShaderRef, 1, &fragSrc, NULL);
	glCompileShader(fragmentShaderRef);

	// Get VERTEX log
	GLint logLength;
	glGetShaderiv(vertexShaderRef, GL_INFO_LOG_LENGTH , &logLength);
	if (logLength > 1)
	{
		printf("Printing log for VERTEX vertex:\n\n");
		GLchar* compiler_log = (GLchar*)malloc(logLength);
		glGetShaderInfoLog(vertexShaderRef, logLength, 0, compiler_log);
		printf("%s\n", compiler_log);
		free (compiler_log);
	}
	// Get FRAG log
	glGetShaderiv(fragmentShaderRef, GL_INFO_LOG_LENGTH , &logLength);
	if (logLength > 1)
	{
		printf("Printing log for FRAGMENT shader:\n\n");
		GLchar* compiler_log = (GLchar*)malloc(logLength);
		glGetShaderInfoLog(fragmentShaderRef, logLength, 0, compiler_log);
		printf("%s\n", compiler_log);
		free (compiler_log);
	}

	// Asserting successful compile
	GLint status;
	glGetShaderiv(vertexShaderRef, GL_COMPILE_STATUS, &status);
	assert(status == GL_TRUE);
	glGetShaderiv(fragmentShaderRef, GL_COMPILE_STATUS, &status);
	assert(status == GL_TRUE);

	// LINKING / MAKING THE PROGRAM
	GLuint shaderProgramRef = glCreateProgram();
	glAttachShader(shaderProgramRef, vertexShaderRef);
	glAttachShader(shaderProgramRef, fragmentShaderRef);
	// (possibly specify which output -> which buffer)
	// └──▶ glBindFragDataLocation(shaderProgram, 0, "outColor");
	glLinkProgram(shaderProgramRef);

	GLint isLinked = 0;
	glGetProgramiv(shaderProgramRef, GL_LINK_STATUS, &isLinked);
	if(isLinked == GL_FALSE)
	{
		printf("Printing log for linking status:\n\n");
		GLint logLength = 0;
		glGetProgramiv(shaderProgramRef, GL_INFO_LOG_LENGTH, &logLength);
		GLchar* linker_log = (GLchar*)malloc(logLength);
		glGetProgramInfoLog(shaderProgramRef, logLength, &logLength, linker_log);
		printf("%s\n", linker_log);
		free (linker_log);
		exit(1);
	}
	std::cout << "Program linked successfully." << std::endl;
	 
	//Always detach shaders after a successful link.
	glDetachShader(shaderProgramRef, vertexShaderRef);
	glDetachShader(shaderProgramRef, fragmentShaderRef);
	glDeleteShader(fragmentShaderRef);
	glDeleteShader(vertexShaderRef);
	return shaderProgramRef;
}

GLuint createVertexArrayObject()
{
	GLuint VAO;
	glGenVertexArrays(1, &VAO);
	return VAO;
}

void setVertexAttribPointer(GLuint shaderProgram, const char* name, GLint numComponents, GLenum type, GLsizei stride, int offset)
{
	GLint attrib = glGetAttribLocation(shaderProgram, name);
	glEnableVertexAttribArray(attrib);
	glVertexAttribPointer(attrib, numComponents, type, GL_FALSE, stride, (void*)((intptr_t)offset));
}


void setSingleFormat(const char *name, int numComponents, char ctype, GLsizei stride, int offset,
					GLuint shaderProgram)
{

	// Read specs
		GLenum type;
		int typeSize;
		switch( ctype ) {
			case 'b': type = GL_BYTE; 			typeSize = sizeof(GLbyte);		break;
			case 'B': type = GL_UNSIGNED_BYTE; 	typeSize = sizeof(GLubyte);		break;
			case 's': type = GL_SHORT;          typeSize = sizeof(GLshort);		break;
			case 'S': type = GL_UNSIGNED_SHORT;	typeSize = sizeof(GLushort);	break;
			case 'i': type = GL_INT;            typeSize = sizeof(GLint);		break;
			case 'I': type = GL_UNSIGNED_INT;   typeSize = sizeof(GLuint);		break;
			case 'f': type = GL_FLOAT;          typeSize = sizeof(GLfloat);		break;
			default:  type = 0;                 assert(!"Unknown type."); break;
		}
	// printf("Setting format..    %s, %d%c, %d, %d", name, numComponents, ctype, stride, offset);
	std::cout << std::endl;

	GLint attrib = glGetAttribLocation(shaderProgram, name);
	glEnableVertexAttribArray(attrib);
	glVertexAttribPointer(attrib, numComponents, type, GL_FALSE, stride * typeSize, (void*)((intptr_t)(offset * typeSize)));
}

/* NOTE: Limit of name = 30 bytes */
struct SingleFormat {
	SingleFormat(int offset, int numComponents, char type, char name[30])
	{
		this->offset = offset; this->numComponents = numComponents; this->ctype = type; strcpy(this->name, name);
	}
	int offset;
	int numComponents;
	char ctype;
	char name[30];
};

void setFormat(const char *format, GLuint shaderProgram)
{
	std::vector<SingleFormat> formats = {};

	char c;
	char *nameBuffer = new char[30];
	int nameIndex = 0;
	
	int totalOffset = 0;

	const char *ptr = format;
	do {
		c = *(ptr); // Store next character in c
		if (c == 0) break;
		if (c == ' ') { // Delimiter - end of name
			// Read next 2 characters
			ptr ++;
			int numComponents = static_cast<int>(*ptr - '0');
			ptr ++;
			char typeChar = *ptr;
			// End name.
			nameBuffer[nameIndex] = 0;
			nameIndex = 0;
			// Create single format object
			SingleFormat newf(totalOffset, numComponents, typeChar, nameBuffer);
			totalOffset += numComponents;
			formats.push_back(newf);
			// Move pointer to next space if any - if 0 is encountered, break whole loop
				bool hasToBreak = false;
				do {
					c = *(++ ptr);
					if (c == 0)  {
						hasToBreak = true;
						break;
					}
				} while (c != ' ');

				if (hasToBreak) break;
			//
		} else { // Add to name:
			nameBuffer[nameIndex] = c;
			nameIndex ++;
		}
		ptr ++;
		
	} while (c != 0);

	// Loop through and handle formats
	for (auto it = formats.begin(); it != formats.end(); it ++)
	{
		setSingleFormat(it->name, it->numComponents, it->ctype, totalOffset, it->offset, shaderProgram);
	}
	
	delete[] nameBuffer;
}



// CALLBACK

void APIENTRY openglCallbackFunction(GLenum source,
                                           GLenum type,
                                           GLuint id,
                                           GLenum severity,
                                           GLsizei length,
                                           const GLchar* message,
                                           const void* userParam)
{
	using namespace std;
    if (type == GL_DEBUG_TYPE_OTHER) return;
 
    cout << "---------------------opengl-callback-start------------" << endl;
    cout << "message: "<< message << endl;
    cout << "type: ";
    switch (type) {
    case GL_DEBUG_TYPE_ERROR:
        cout << "ERROR";
        break;
    case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
        cout << "DEPRECATED_BEHAVIOR";
        break;
    case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
        cout << "UNDEFINED_BEHAVIOR";
        break;
    case GL_DEBUG_TYPE_PORTABILITY:
        cout << "PORTABILITY";
        break;
    case GL_DEBUG_TYPE_PERFORMANCE:
        cout << "PERFORMANCE";
        break;
    case GL_DEBUG_TYPE_OTHER:
        cout << "OTHER";
        break;
    }
    cout << endl;
 
    cout << "id: " << id << endl;
    cout << "severity: ";
    switch (severity){
    case GL_DEBUG_SEVERITY_LOW:
        cout << "LOW";
        break;
    case GL_DEBUG_SEVERITY_MEDIUM:
        cout << "MEDIUM";
        break;
    case GL_DEBUG_SEVERITY_HIGH:
        cout << "HIGH";
        break;
    }
    cout << endl;
    cout << "---------------------opengl-callback-end--------------" << endl;
}




// BOILER PLATE
void GLFW_boilerPlate(GLFWwindow **window, GLFWerrorfun error_callback)
{
	glfwSetErrorCallback(error_callback);
	if (!glfwInit())
		exit(EXIT_FAILURE);
	// Set hints
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	// For debugging
	glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
	

	// Create window and context
	*window = glfwCreateWindow(640, 480, "Program", NULL, NULL);
	if (!window) {
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(*window);
	// Initialize GLEW
	glewExperimental = GL_TRUE;
	if (glewInit()) fprintf(stderr, "Failed to initialize GLEW.\n");

	glDebugMessageCallback(openglCallbackFunction, 0);
}

