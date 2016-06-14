CC = g++
CFLAGS = -c -std=c++11 -O3
LDFLAGS =  -L/usr/local/lib/

CUDA_PATH = /usr/local/cuda/

COMPUTE_CAPABILITY := \
-gencode arch=compute_20,code=sm_20 \
-gencode arch=compute_30,code=sm_30 \
-gencode arch=compute_50,code=sm_50 \
-gencode arch=compute_52,code=sm_52

NVCC = $(CUDA_PATH)/bin/nvcc
NVCCFLAGS = $(CFLAGS) $(COMPUTE_CAPABILITY) -D_FORCE_INLINES -DGPU

# Makefile for gerbil project
GPU = true

# Source Files
CU_SRCS += \
src/cuda_ds/AddKernel.cu \
src/cuda_ds/CompressKernel.cu

CPP_SRCS += \
src/gerbil/Application.cpp \
src/gerbil/Bundle.cpp \
src/gerbil/FastFile.cpp \
src/gerbil/FastParser.cpp \
src/gerbil/FastReader.cpp \
src/gerbil/KmcWriter.cpp \
src/gerbil/KmerDistributor.cpp \
src/gerbil/SequenceSplitter.cpp \
src/gerbil/SuperReader.cpp \
src/gerbil/SuperWriter.cpp \
src/gerbil/TempFile.cpp \
src/gerbil/debug.cpp \
src/gerbil/gerbil.cpp

CPP_OBJECTS := $(CPP_SRCS:.cpp=.o) 
CU_OBJECTS := $(CU_SRCS:.cu=.o)

SRCS = $(CPP_SRCS)
OBJECTS = $(CPP_OBJECTS)

ifeq ($(GPU),true)
	CFLAGS += -DGPU -I$(CUDA_PATH)/include/
	SRCS += $(CU_SRCS)
	OBJECTS += $(CU_OBJECTS)
	LDFLAGS = -L$(CUDA_PATH)/lib64 -lcudart
endif

LIBS := -lboost_system -lboost_thread -lboost_filesystem -lboost_regex -lbz2 -lz -lpthread

RM := rm -rf
EXECUTABLE := gerbil

all: $(SRCS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o  $(EXECUTABLE) $(OBJECTS) $(LDFLAGS) $(LIBS)
	
%.o: %.cu
	$(NVCC) $(CFLAGS) -o "$@" "$<"

%.o: %.cpp
	$(CC) $(CFLAGS) -o "$@" "$<"


# Other Targets
clean:
	-$(RM) $(OBJECTS) $(CU_OBJECTS) $(EXECUTABLE) prototype

.PHONY: all clean
