INCLUDE_DIRECTORIES(/data1/lz/LilianTests/ARPACKLAB/arpack++/include/)
INCLUDE_DIRECTORIES(/data1/lz/LilianTests/ARPACKLAB/arpack++)

#project name
PROJECT(LILIANTESTS)

cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET(BUILD_SHARED_LIBS ON)

## where are user-specific cmake modules
SET(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})

OPTION(USE_BIAS OFF)

IF(USE_BIAS)
FIND_PACKAGE(BIAS)
IF(BIAS_FOUND)
  INCLUDE(${BIAS_USE_FILE})
ELSE(BIAS_FOUND)
  MESSAGE(SEND_ERROR "BIAS lib not found.")
ENDIF(BIAS_FOUND)
ENDIF(USE_BIAS)


IF(USE_BIAS)
INCLUDE_DIRECTORIES(${BIAS_INCLUDE_DIR} ${BIAS_INCLUDE_DIR}/BIAS ${WXWIDGETS_INCLUDE_DIR})
LINK_DIRECTORIES( ${BIAS_LINK_DIRECTORIES})
ENDIF(USE_BIAS)
  
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BIAS_CXX_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  ${BIAS_C_FLAGS}")

# source files of library "LineMatchingLib" to be created
SET(LineMatchingLib_SRCS
    PairwiseLineMatching.cpp
    LineDescriptor.cpp
    EDLineDetector.cpp
)
# header files to be installed
SET(LineMatchingLib_HEADER
    PairwiseLineMatching.hh
    LineDescriptor.hh
    EDLineDetector.hh
    LineStructure.hh
)


ADD_LIBRARY(LineMatchingLib
            ${LineMatchingLib_SRCS}
            ${LineMatchingLib_HEADER})
TARGET_LINK_LIBRARIES(LineMatchingLib  arpack /data1/lz/LilianTests/ARPACKLAB/SuperLU/libsuperlu.so  ${BIAS_LIBRARIES}) 


ADD_EXECUTABLE(TestLineMatchingAlgorithm TestLineMatchingAlgorithm.cpp)
TARGET_LINK_LIBRARIES(TestLineMatchingAlgorithm   LineMatchingLib)


