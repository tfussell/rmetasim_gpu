
# User flags
WITHGPU = 0
DEBUG = 1
VERBOSE = 1
ARCH := sm_21


# Locations
GPUSRCDIR = gpu


# Enable/disable verbosity by preprending @
ifeq ($(VERBOSE), 1)
VPREFIX := 
else
VPREFIX := @
endif


# Compilers/Linkers
NVCC := $(VPREFIX)/usr/local/cuda/bin/nvcc
CXX  := $(VPREFIX)/usr/bin/g++-4.4
RM   := $(VPREFIX)rm

ifeq ($(WITHGPU),1)
LINKER := $(NVCC)
else
LINKER := $(CXX)
endif


# Warning flags
NVCCWARNFLAGS := 
CXXWARNFLAGS  := -W -Wall -Wextra -Wvla -Wuninitialized -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -pedantic -Werror


# External Dependencies
INCLUDES  := . /usr/local/lib/R/include /usr/local/include
LIBRARIES := gfortran quadmath m Rblas Rlapack R
LIBDIRS   := /usr/local/lib /usr/local/lib/R/lib /usr/X11R6/lib64 /usr/lib64 /home/thomas/Download/R-2.14.1/src/main

ifeq ($(WITHGPU),1)
INCLUDES  += /usr/local/cuda/include 
LIBRARIES += cuda GL GLU X11 Xi Xmu glut cufft cublas curand cusparse 
LIBDIRS   += /usr/lib/nvidia-current
endif

# Compiler/Linker Flags
NVCCFLAGS := --compiler-options -fPIC -D__Linux__ -arch=$(ARCH) $(addprefix -I, $(INCLUDES)) $(NVCCWARNFLAGS)
CXXFLAGS  := -fPIC -D__Linux__ $(addprefix -isystem, $(INCLUDES)) $(CXXWARNFLAGS)
LINKFLAGS := -shared $(addprefix -l, $(LIBRARIES)) $(addprefix -L, $(LIBDIRS))


# Debug Options
ifeq ($(DEBUG),1)
NVCCFLAGS += -D_DEBUG -g
CXXFLAGS  += -D_DEBUG -ggdb
else 
NVCCFLAGS += --compiler-options -fno-strict-aliasing --compiler-options -use_fast_math
CXXFLAGS  += -fno-strict-aliasing
endif


# Objects to be compiled
OBJECTS    := AlleleObj.o AlleleTbl.o BaseObj.o Democlass.o Eigen.o FastAllele.o FastSeqAllele.o Genotype.o Landscape.o PackedIndividual.o rmetasim.o RandLib.o SiteObj.o TransMat.o
GPUOBJECTS := CudaKernels.o Landscape_gpu.o PRNG.o rmetasim_gpu.o SimulationTimer.o

ifeq ($(WITHGPU),1)
OBJECTS += $(GPUOBJECTS)
endif


# Final Binary
TARGET = rmetasim.so


# Canned compilation commands
define cpp-compile
$(CXX) $(CXXFLAGS) -o $@ -c $<
endef

define nvcc-compile
$(NVCC) $(NVCCFLAGS) -o $@ -c $<
endef


# Rules
all: $(TARGET) Makefile

$(TARGET): $(OBJECTS) Makefile
	@echo "Linking objects into" $(TARGET)	
	$(LINKER) $(LINKFLAGS) $(OBJECTS) -o $(TARGET)

%.o : %.cc Makefile
	@echo "Compiling:" $@
	$(cpp-compile)

ifeq ($(WITHGPU),1)
%.o : $(GPUSRCDIR)/%.cu Makefile
	@echo "Compiling:" $@
	$(nvcc-compile)

%.o : $(GPUSRCDIR)/%.cc Makefile
	@echo "Compiling:" $@
	$(cpp-compile)
endif

clean:
	@echo "Cleaning all objects"
	$(RM) -f $(OBJECTS) $(TARGET)

.PHONY: all clean
