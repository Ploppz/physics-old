SHELL = /bin/bash
CXX = clang++
OBJDIR = .obj
TOOLDIR = ~/scripts/make
DEFAULT_MAIN = src/main.c++
# We have to add the include directory because freetype uses relative (to freetype2 dir) include file names
FLAGS = -std=c++11 -ggdb -iquote src -I/usr/include/freetype2 -Wall -Wno-sign-compare -Wno-unused-variable 
# FLAGS = -std=c++11 -ggdb -Isrc -I/usr/include/freetype2 -Wall -Wextra
LDFLAGS = -lglfw -lGL -lGLEW -lfreetype
DEFINE = -DFREETYPE_GL_USE_VAO -DUSE_GLFW

MAIN := $(DEFAULT_MAIN)
SOURCE_CPP := $(shell $(TOOLDIR)/find-filetype src c++) $(MAIN)
OBJECTS := $(shell $(TOOLDIR)/getobjects src c++ c)

# Add object file of the main file
# OBJECTS += $(patsubst */%.c++, .obj/%.o, $(MAIN))
MAIN_STEM := $(shell echo $(MAIN) | sed 's/.*\/\(.*\)\.c++/\1/g')
OBJECTS += .obj/$(MAIN_STEM).o
ifneq "$(MAIN)" "$(DEFAULT_MAIN)"
ifeq "$(USE_CATCH)" "true"
	OBJECTS += .obj/catch.o
endif
endif

PROGRAM = bin/$(MAIN_STEM)

SHADERS := $(shell $(TOOLDIR)/find-filetype src glsl)
FRAGSHADERS := $(shell $(TOOLDIR)/find-filetype src frag)
VERTSHADERS := $(shell $(TOOLDIR)/find-filetype src vert)
SHADER_OBJ := $(shell $(TOOLDIR)/getshaderobjects src)
MSG = $(TOOLDIR)/msg

#--------------------------------------------

COMMA=","

# TODO progress:
# .depend doesn't work explicitly like wanted...

$(PROGRAM) : objects

objects : $(OBJECTS) $(SHADER_OBJ)
	@$(MSG) "Link"
	$(CXX) -o $(PROGRAM) $(OBJECTS) $(SHADER_OBJ) $(LDFLAGS) $(FLAGS)
	#
	rm -f $(patsubst src/%.glsl, src/%.c++, $(SHADERS))
	@echo -e '\e[0m'

# Source files
$(OBJDIR)/%.o : src/%.c++ | $(OBJDIR) .depend
	@$(MSG) "Compile $<"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS)
	@mkdir -p $$(dirname "$@")

$(OBJDIR)/%.o : src/%.c | $(OBJDIR) .depend
	@$(MSG) "Compile $<"
	$(CXX) -c $< -o $@ $(DEFINE)

# Unit tests main files
$(OBJDIR)/%.o : tests/%.c++ | $(OBJDIR) .depend
	@$(MSG) "Compile unit test main file - $<"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS)
	@mkdir -p $$(dirname "$@")

# SHADER GENERATION
$(OBJDIR)/%.o : src/%.glsl | $(OBJDIR) .depend
	@$(MSG) "$*.glsl ---> $*.c++ ---> $*.o"
	@echo 'namespace shaders { extern char const * const $(shell $(TOOLDIR)/deslash $*) = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*.c++
	@$(CXX) -c src/$*.c++ -o $@ $(DEFINE) $(FLAGS)
	@rm src/$*.c++

$(OBJDIR)/%_v.o : src/%.vert | $(OBJDIR) .depend
	@$(MSG) "$*.vert ---> $*.c++ ---> $*.o"
	@echo 'namespace shaders { extern char const * const $(shell $(TOOLDIR)/deslash $*)_v = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*_v.c++
	@$(CXX) -c src/$*_v.c++ -o $@ $(DEFINE) $(FLAGS)
	@rm src/$*_v.c++

$(OBJDIR)/%_f.o : src/%.frag | $(OBJDIR) .depend
	@$(MSG) "$*.frag ---> $*.c++ ---> $*.o"
	@echo 'namespace shaders { extern char const * const $(shell $(TOOLDIR)/deslash $*)_f = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*_f.c++
	@$(CXX) -c src/$*_f.c++ -o $@ $(DEFINE) $(FLAGS)
	@rm src/$*_f.c++


# src/shaders.h : $(SHADERS) $(FRAGSHADERS) $(VERTSHADERS)
src/shaders.h :
	@$(MSG) "Generate src/shaders.h"
	@echo "namespace shaders {" > src/shaders.h
	@$(TOOLDIR)/getshadernames src | sed 's/\(.*\)/\textern char const * const \1;/g' >> src/shaders.h
	@echo "}" >> src/shaders.h

# Order-only targets
$(OBJDIR):
	@$(MSG) "Create directory structure of $(OBJDIR)"
	@mkdir -p $(OBJDIR)
	@cd $(OBJDIR); (cd ../src; find -type d ! -name .) | xargs mkdir
	@mkdir $(OBJDIR)/tests

.depend:
	@$(MSG) "Generate dependencies"
	$(TOOLDIR)/gendepend $(SOURCE_CPP) > .depend


# Phony targets
.PHONY: clean
clean:
	rm -rf $(PROGRAM) $(OBJDIR) $(patsubst src/%.glsl, src/%.c++, $(SHADERS)) $(patsubst src/%.frag, src/%.c++, $(FRAGSHADERS)) $(patsubst src/%.vert, src/%.c++, $(VERTSHADERS)) src/shaders.h .depend

.PHONY: flags
flags:
	@echo "$(FLAGS)"

-include .depend
