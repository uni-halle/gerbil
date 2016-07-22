cmake_minimum_required(VERSION 2.8)
project(Gerbil)

# The version number.
set (Gerbil_VERSION_MAJOR 0)
set (Gerbil_VERSION_MINOR 2)

# Compile Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

# Include Files
set(HEADER_FILES
        include/gerbil/Application.h
        include/gerbil/Bundle.h
        include/gerbil/CpuHasher.h
        include/gerbil/debug.h
        include/gerbil/FailureBuffer.h
        include/gerbil/FastFile.h
        include/gerbil/FastParser.h
        include/gerbil/FastReader.h
        include/gerbil/GpuHasher.h
        include/gerbil/KmcWriter.h
        include/gerbil/KMer.h
        include/gerbil/KmerCountingHashTable.h
        include/gerbil/KmerDistributer.h
        include/gerbil/KmerHasher.h
        include/gerbil/SequenceSplitter.h
        include/gerbil/SuperReader.h
        include/gerbil/SuperWriter.h
        include/gerbil/SyncQueue.h
        include/gerbil/TempFile.h
        include/gerbil/TempFileStatistic.h
        include/gerbil/ThreadBarrier.h
        include/gerbil/types.h
        include/cuda_ds/CountingHashTable.h
        include/cuda_ds/CountingHashTable/_CountingHashTable.h
        include/cuda_ds/CountingHashTable/_Key.h
        include/cuda_ds/CountingHashTable/_KeyValuePair.h
        )

# Source Files
set(SOURCE_FILES_CPP
        src/gerbil/Application.cpp
        src/gerbil/Bundle.cpp
        src/gerbil/FastFile.cpp
        src/gerbil/FastParser.cpp
        src/gerbil/FastReader.cpp
        src/gerbil/KmcWriter.cpp
        src/gerbil/KmerDistributor.cpp
        src/gerbil/SequenceSplitter.cpp
        src/gerbil/SuperReader.cpp
        src/gerbil/SuperWriter.cpp
        src/gerbil/TempFile.cpp
        src/gerbil/debug.cpp
        src/gerbil/gerbil.cpp
        )

set(SOURCE_FILES_CUDA
        src/cuda_ds/AddKernel.cu
        src/cuda_ds/CompressKernel.cu
        )

# import mandatory libraries
find_package (Boost COMPONENTS REQUIRED system thread filesystem regex)
find_package (BZip2 REQUIRED)
find_package (ZLIB REQUIRED)

set(MANDATORY_INCLUDES ${Boost_INCLUDE_DIRS} ${BZip2_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
set(MANDATORY_LIBS  ${Boost_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES})

# import optional libraries
find_package(CUDA)
if(CUDA_FOUND)
    set(GPU 1)                                              # set GPU flag to true in config file
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)                      # turn off flag propagation
    set(EXTRA_INCLUDES ${CUDA_INCLUDE_DIRS})                # use cuda headers
    set (EXTRA_LIBS ${EXTRA_LIBS} ${CUDA_LIBRARIES})        # use cuda libraries
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}                  # set nvcc compile flags
            -D_FORCE_INLINES -std=c++11 -O3)
    cuda_compile(cuda_objs ${SOURCE_FILES_CUDA})            # compile .cu sources
endif(CUDA_FOUND)

# configure a header file to pass some of the CMake settings to the source code
configure_file (
        "${PROJECT_SOURCE_DIR}/include/gerbil/config.h.in"
        "${PROJECT_SOURCE_DIR}/include/gerbil/config.h"
)

# include headers
include_directories(${MANDATORY_INCLUDES} ${EXTRA_INCLUDES})

# build executable
add_executable(Gerbil ${HEADER_FILES} ${SOURCE_FILES_CPP} ${cuda_objs})

# link against external libraries
target_link_libraries (Gerbil ${MANDATORY_LIBS} ${EXTRA_LIBS})

# install binary
install (TARGETS Gerbil DESTINATION bin)

