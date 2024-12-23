cmake_minimum_required(VERSION 3.10)
project(FluidSimulatorProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT DEFINED TYPES)
    set(TYPES "FLOAT,FIXED(32,16),DOUBLE")
endif()

if(NOT DEFINED SIZES)
    set(SIZES "S(10,10),S(1920,1080),S(36,84)")
endif()

add_definitions(-DTYPES="${TYPES}" -DSIZES="${SIZES}")

include_directories(include)

set(SOURCES
    src/main.cpp
    src/utils.cpp
)

add_executable(FluidSimulatorExecutable ${SOURCES})