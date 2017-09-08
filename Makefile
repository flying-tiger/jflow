#-----------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------
VPATH    := source test
CXX      := cl.exe
CXXFLAGS := /I./extern/catch/single_include \
            /experimental:module            \
            /std:c++latest                  \
            /module:search obj/             \
            /nologo                         \
            /EHsc                           \
            /MD


#-----------------------------------------------------------------------
# Word Lists
#-----------------------------------------------------------------------
TEST_FILES := $(notdir $(wildcard test/*.cpp))
TEST_OBJ   := $(addprefix obj/, $(TEST_FILES:.cpp=.obj))

#-----------------------------------------------------------------------
# Targets
#-----------------------------------------------------------------------

.PHONY: default
default: jflow

.PHONY: jflow
jflow: bin/jflow.exe

.PHONY: test
test: bin/test_jflow.exe
	./$<

.PHONY: clean
clean:
	rm -rf obj bin



#-----------------------------------------------------------------------
# Dependency Graphs
#-----------------------------------------------------------------------

bin/jflow.exe: obj/main.obj
obj/main.obj:  obj/jflow.euler.ifc

bin/test_jflow.exe: $(TEST_OBJ)
obj/test_euler.obj: obj/jflow.euler.ifc


#-----------------------------------------------------------------------
# Build Rules
#-----------------------------------------------------------------------
bin/%: | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fe./$@ $^

obj/%.obj: %.cpp | obj
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo./obj/ /c $<

obj/jflow.%.ifc: %.ixx | obj
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo./obj/ /module:output $@ /c $<

obj bin:
	mkdir -p $@
