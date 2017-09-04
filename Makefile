#-----------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------
VPATH    := source test
CXX      := cl.exe
CXXFLAGS := /I./extern/catch/single_include \
            /experimental:module            \
            /std:c++latest                  \
            /nologo                         \
            /EHsc                           \
            /MD                             \
            /W4 /WX


#-----------------------------------------------------------------------
# Word Lists
#-----------------------------------------------------------------------
MODULES := euler


#-----------------------------------------------------------------------
# Targets
#-----------------------------------------------------------------------
.PHONY: clean
clean:
	rm -rf obj bin

.PHONY: modules
modules: $(addprefix obj/, $(addsuffix .ifc, $(MODULES)))

.PHONY: test
test: bin/test.exe
	./$<

.PHONY: check
check:
	@echo $(TEST_FILES)
	@echo $(TEST_OBJ)


#-----------------------------------------------------------------------
# Build Rules
#-----------------------------------------------------------------------

TEST_FILES := $(notdir $(wildcard test/*.cpp))
TEST_OBJ   := $(addprefix obj/, $(TEST_FILES:.cpp=.obj))
bin/test.exe: $(TEST_OBJ)

bin/%: | bin
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fe./$@ $^

obj/%.obj: %.cpp | obj
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo./obj/ /c $<

obj/%.ifc: %.ixx | obj
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) /Fo./obj/ /module:output $@ -c $<

obj bin:
	mkdir -p $@
