SHELL := /bin/bash
# OBJECTS = $(patsubst src/%.cpp, .obj/%.o, $(wildcard src/*.cpp)) $(patsubst src/%.c, .obj/%.o, $(wildcard src/*.c))
MSG = ~/scripts/make/msg
OBJDIR = .obj
SOURCE_CPP = $(shell ~/scripts/make/find-filetype src cpp)
OBJECTS = $(shell /home/ploppz/scripts/make/getobjects src cpp c)

SHADERS = $(shell ~/scripts/make/find-filetype src glsl)
FRAGSHADERS = $(shell ~/scripts/make/find-filetype src frag)
VERTSHADERS = $(shell ~/scripts/make/find-filetype src vert)
SHADER_OBJ = $(shell ~/scripts/make/getshaderobjects src)
PROGRAM = bin/program
DEFINE = -DFREETYPE_GL_USE_VAO -DUSE_GLFW
# We have to add the include directory because freetype uses relative (to freetype2 dir) include file names
FLAGS = -ggdb -I/usr/include/freetype2
LDFLAGS = -lglfw -lGL -lGLEW -lfreetype
CXX = g++

#--------------------------------------------

COMMA=","

$(PROGRAM) : objects

objects : $(OBJECTS) $(SHADER_OBJ)
	@$(MSG) "Link"
	$(CXX) -o $(PROGRAM) $(OBJECTS) $(SHADER_OBJ) $(LDFLAGS) -std=c++11
	#
	rm -f $(patsubst src/%.glsl, src/%.cpp, $(SHADERS))
	#

$(OBJDIR)/%.o : src/%.cpp | $(OBJDIR)
	@$(MSG) "cpp - $<"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS) -std=c++11
	@mkdir -p $$(dirname "$@")

$(OBJDIR)/%.o : src/%.c | $(OBJDIR)
	@$(MSG) "c"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS) -std=c++11

# SHADER GENERATION
$(OBJDIR)/%.o : src/%.glsl | $(OBJDIR)
	@$(MSG) "GENERATE & COMPILE SHADER CPP $*.cpp	from	$*.glsl"
	@echo 'namespace shaders { extern char const * const $(shell ~/scripts/make/deslash $*) = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*.cpp
	@$(CXX) -c src/$*.cpp -o $@ $(DEFINE) $(FLAGS) -std=c++11
	@rm src/$*.cpp

$(OBJDIR)/%_v.o : src/%.vert | $(OBJDIR)
	@echo 'namespace shaders { extern char const * const $(shell ~/scripts/make/deslash $*)_v = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*_v.cpp
	@$(CXX) -c src/$*_v.cpp -o $@ $(DEFINE) $(FLAGS) -std=c++11
	@rm src/$*_v.cpp

$(OBJDIR)/%_f.o : src/%.frag | $(OBJDIR)
	@echo 'namespace shaders { extern char const * const $(shell ~/scripts/make/deslash $*)_f = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*_f.cpp
	@$(CXX) -c src/$*_f.cpp -o $@ $(DEFINE) $(FLAGS) -std=c++11
	@rm src/$*_f.cpp


src/shaders.h : $(SHADERS) $(FRAGSHADERS) $(VERTSHADERS)
	##GENERATE SHADER HEADER shaders.h
	@echo "namespace shaders {" > src/shaders.h
	@~/scripts/make/getshadernames src | sed 's/\(.*\)/\textern char const * const \1;/g' >> src/shaders.h
	@echo "}" >> src/shaders.h

$(OBJDIR):
	@$(MSG) "Create directory structure of $(OBJDIR)"
	@mkdir -p $(OBJDIR)
	@cd $(OBJDIR); (cd ../src; find -type d ! -name .) | xargs mkdir


# note: ensures that .depend is unconditionally updated??
depend : .depend

.depend : $(SOURCE_CPP)
	@$(MSG) ".depend"
	@# DID: don't cd, then you don't need to prefix with src (?)
	@ # Prefix prereqs with src/, and objects with $(OBJDIR)/
	$(CXX) -MM -MG $^ -std=c++11 | sed 's/\(.*\:\)/$(OBJDIR)\/\1/g' > .depend;

clean:
	rm -rf $(PROGRAM) $(OBJDIR) $(patsubst src/%.glsl, src/%.cpp, $(SHADERS)) $(patsubst src/%.frag, src/%.cpp, $(FRAGSHADERS)) $(patsubst src/%.vert, src/%.cpp, $(VERTSHADERS)) src/shaders.h .depend

-include .depend
