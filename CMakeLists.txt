# CMake main for simul++
cmake_minimum_required(VERSION 3.1...3.22)

# Compiler specification (not needed)
set (CMAKE_CXX_COMPILER "g++")

# Project name and version 
project (
    simul++
    VERSION 1.5
    LANGUAGES CXX
)

string(TIMESTAMP DATE "%d.%m.%Y %H:%M")

configure_file (
    "${simul++_SOURCE_DIR}/main/incl/version.hpp.in"
    "${simul++_SOURCE_DIR}/main/incl/version.hpp"
)

option(MYDEBUG "Sets CMAKE_BUILD_TYPE to DEBUG and adds #define DEBUG and -Wall" off)
option(NOPROGRESSBAR "Sets NOPROGRESSBAR to defined, no progress bar is printed to stdout..." off)
option(VERLET1 "Velocity Verlet version of velocity for Verlet-type integrators" on)
option(VERLET2 "Harmonic average for velocity (v^2 = v(i-1/2)*v(i+1/2)) for Verlet-type integrators" off)
option(VERLET0 "Velocity v(i+1/2) for Verlet-type integrators" off)
option(VERLET3 "Velocity as a mean of both shifted energies (E(i-1/2), E(i+1/2)) for Verlet-type integrators" off)
option(FASTMATH "Enables -ffast-math optimization flag (cca 10 % speedup), can be unsafe; not tested" off)
option(PARALLEL "Parallelization using OpenMP" off)

set( CMAKE_VERBOSE_MAKEFILE off )

# Output directory for executables
if (MYDEBUG)
set (CMAKE_BUILD_TYPE DEBUG)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${simul++_SOURCE_DIR}/bin-debug)
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -Wextra -O3 -g")
else ()
set (CMAKE_BUILD_TYPE RELEASE)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${simul++_SOURCE_DIR}/bin)
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3")
endif ()

# Progress bar option
if (NOPROGRESSBAR)
message ("simul++ will be compiled without progress bar")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOPROGRESSBAR")
endif()

# Velocity version for Verlet-type integrators
if (VERLET0)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERLET=0")
elseif (VERLET2)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERLET=2")
elseif (VERLET3)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERLET=3")
elseif (VERLET1)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERLET=1")
else ()
message ("VERLET variable not set, using default VERLET=1")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERLET=1")
endif ()

# Enable -ffast-math
if (FASTMATH)
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
message ("Warning: -ffast-math used, can be unsafe!")
endif ()

# Parallel - use OpenMP
if (PARALLEL)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARALLEL -fopenmp")
endif ()

# Print options summary
message("BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")

# Add subdirectories for libraries and main
add_subdirectory (main)
add_subdirectory (mymath)
add_subdirectory (system)
add_subdirectory (utils)
add_subdirectory (integrators)
