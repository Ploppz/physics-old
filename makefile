SHELL = /bin/bash
DEFAULT_MAIN = src/main.c++
MAIN := $(DEFAULT_MAIN)
MSG = tools/msg
OBJDIR = .obj
SOURCE_CPP := $(shell ~/scripts/make/find-filetype src c++) $(MAIN)
OBJECTS := $(shell /home/ploppz/scripts/make/getobjects src c++ c)

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

SHADERS := $(shell ~/scripts/make/find-filetype src glsl)
FRAGSHADERS := $(shell ~/scripts/make/find-filetype src frag)
VERTSHADERS := $(shell ~/scripts/make/find-filetype src vert)
SHADER_OBJ := $(shell ~/scripts/make/getshaderobjects src)
# We have to add the include directory because freetype uses relative (to freetype2 dir) include file names
FLAGS = -std=c++11 -ggdb -Isrc -I/usr/include/freetype2 -Wall -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable
# FLAGS = -std=c++11 -ggdb -Isrc -I/usr/include/freetype2 -Wall -Wextra
LDFLAGS = -lglfw -lGL -lGLEW -lfreetype
DEFINE = -DFREETYPE_GL_USE_VAO -DUSE_GLFW
CXX = g++

#--------------------------------------------

COMMA=","

# TODO progress:
# Fixed the shaders.h problem BUT files compile before Make generates src/shaders.h

$(PROGRAM) : objects

objects : $(OBJECTS) $(SHADER_OBJ)
	@$(MSG) "Link"
	$(CXX) -o $(PROGRAM) $(OBJECTS) $(SHADER_OBJ) $(LDFLAGS) $(FLAGS)
	#
	rm -f $(patsubst src/%.glsl, src/%.c++, $(SHADERS))
	@echo -e '\e[0m'

# Source files
$(OBJDIR)/%.o : src/%.c++ | $(OBJDIR)
	@$(MSG) "Compile $<"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS)
	@mkdir -p $$(dirname "$@")

$(OBJDIR)/%.o : src/%.c | $(OBJDIR)
	@$(MSG) "Compile $<"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS)

# Unit tests main files
$(OBJDIR)/%.o : tests/%.c++ | $(OBJDIR)
	@$(MSG) "Compile unit test main file - $<"
	$(CXX) -c $< -o $@ $(DEFINE) $(FLAGS)
	@mkdir -p $$(dirname "$@")

# SHADER GENERATION
$(OBJDIR)/%.o : src/%.glsl | $(OBJDIR)
	@$(MSG) "GENERATE & COMPILE SHADER CPP $*.c++	from	$*.glsl"
	@echo 'namespace shaders { extern char const * const $(shell ~/scripts/make/deslash $*) = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*.c++
	@$(CXX) -c src/$*.c++ -o $@ $(DEFINE) $(FLAGS)
	@rm src/$*.c++

$(OBJDIR)/%_v.o : src/%.vert | $(OBJDIR)
	@echo 'namespace shaders { extern char const * const $(shell ~/scripts/make/deslash $*)_v = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*_v.c++
	@$(CXX) -c src/$*_v.c++ -o $@ $(DEFINE) $(FLAGS)
	@rm src/$*_v.c++

$(OBJDIR)/%_f.o : src/%.frag | $(OBJDIR)
	@echo 'namespace shaders { extern char const * const $(shell ~/scripts/make/deslash $*)_f = R"SHADER_D('"$$(cat $<)"')SHADER_D";}' > src/$*_f.c++
	@$(CXX) -c src/$*_f.c++ -o $@ $(DEFINE) $(FLAGS)
	@rm src/$*_f.c++


src/shaders.h : $(SHADERS) $(FRAGSHADERS) $(VERTSHADERS)
	##GENERATE SHADER HEADER shaders.h
	@echo "namespace shaders {" > src/shaders.h
	@~/scripts/make/getshadernames src | sed 's/\(.*\)/\textern char const * const \1;/g' >> src/shaders.h
	@echo "}" >> src/shaders.h

$(OBJDIR):
	@$(MSG) "Create directory structure of $(OBJDIR)"
	@mkdir -p $(OBJDIR)
	@cd $(OBJDIR); (cd ../src; find -type d ! -name .) | xargs mkdir
	@mkdir $(OBJDIR)/tests

.depend : $(SOURCE_CPP)
	@# Needs to be silenced so that `flags` can run 'undisturbed'
	@ # Prefix prereqs with src/, and objects with $(OBJDIR)/.
	@# Make sure shaders.h is prefixed with src/shaders.h
	@$(CXX) -MM -MG $^ -std=c++11 | sed 's/\(.*\:\)/$(OBJDIR)\/\1/g' | sed 's/ \(shaders\.h\|\.\.\/shaders\.h\) / src\/shaders.h /g' > .depend;


clean:
	rm -rf $(PROGRAM) $(OBJDIR) $(patsubst src/%.glsl, src/%.c++, $(SHADERS)) $(patsubst src/%.frag, src/%.c++, $(FRAGSHADERS)) $(patsubst src/%.vert, src/%.c++, $(VERTSHADERS)) src/shaders.h .depend
flags:
	@echo "$(FLAGS)"

-include .depend
