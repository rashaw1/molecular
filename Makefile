# Project name
PROJECT = molecular
# Compiler
CXX = icpc

# Run Options
COMMANDLINE_OPTIONS = 

# MKL libraries
MKLROOT = /opt/intel/mkl

# Compiler options
DEBUG = -g -Wall -O0 -std=c++11 -D_GLIBCXX_DEBUG 
OPTIM = -O3 -Wall -std=c++11 -qopenmp -DMKL_LP64 -mkl -DEIGEN_USE_MKL_ALL
COMPILE_OPTIONS = $(DEBUG)

# Header include directories
HEADERS = -I./inc -I/usr/local/Cellar/boost/1.61.0_1/include -I/usr/local/Cellar/eigen/3.2.9/include/eigen3 -I${MKLROOT}/include

# Libraries for linking
LIBS =  -L/usr/local/Cellar/boost/1.61.0_1/lib -lboost_system -lboost_timer  -std=c++11 -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -qopenmp -I${MKLROOT}/include

# Dependency options
DEPENDENCY_OPTIONS = -MM

#-- Do not edit below this line --

# Subdirs to search for additional source files
SUBDIRS := $(shell ls -F | grep "\/")
DIRS := ./ $(SUBDIRS)
SOURCE_FILES := $(foreach d, $(DIRS), $(wildcard $(d)*.cpp) )

# Create an object for every cpp file
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Dependencies
DEPENDENCIES := $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Create .d files
%.d: %.cpp
	$(CXX) $(DEPENDENCY_OPTIONS) $< -MT "$*.o $*.d" -MF $*.d

# Make $(PROJECT) the default target
all: $(DEPENDENCIES) $(PROJECT)

$(PROJECT): $(OBJECTS)
	$(CXX) -o $(PROJECT) $(OBJECTS) $(LIBS)

# Include dependencies

# Compile every cpp file to an object
%.o: %.cpp
	$(CXX) -c $(COMPILE_OPTIONS) -o $@ $< $(HEADERS)

# Build and run project
run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

# Clean and debug
.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(PROJECT) $(DEPENDENCIES)

clean-all: clean depclean
