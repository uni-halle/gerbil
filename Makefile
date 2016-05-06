CC = g++
NVCC = /usr/local/cuda/bin/nvcc
LDFLAGS =  -L/usr/local/lib/

COMPUTE_CAPABILITY := \
-gencode arch=compute_20,code=sm_20 \
-gencode arch=compute_30,code=sm_30 \
-gencode arch=compute_50,code=sm_50 \
-gencode arch=compute_52,code=sm_52

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

CFLAGS = -c -std=c++11 -O3

SRCS = $(CPP_SRCS)
OBJECTS = $(CPP_OBJECTS)

ifeq ($(GPU),true)
	CC = $(NVCC)
	CFLAGS += -DGPU $(COMPUTE_CAPABILITY) -D_FORCE_INLINES
	SRCS += $(CU_SRCS)
	OBJECTS += $(CU_OBJECTS) 
endif

LIBS := -lboost_system -lboost_thread -lboost_filesystem -lboost_regex -lbz2 -lz -lpthread

RM := rm -rf
EXECUTABLE := gerbil

all: $(SRCS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o  $(EXECUTABLE) $(OBJECTS) $(LDFLAGS) $(LIBS)
	
%.o: %.cu
	$(CC) $(CFLAGS) -o "$@" "$<"

%.o: %.cpp
	$(CC) $(CFLAGS) -o "$@" "$<"


# Other Targets
clean:
	-$(RM) $(OBJECTS) $(EXECUTABLE) prototype

.PHONY: all clean
