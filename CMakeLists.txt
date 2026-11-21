cmake_minimum_required(VERSION 3.21)
project(ustsolver VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type configuration
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Build directory: ${PROJECT_BINARY_DIR}")

# User options
option(UST_ENABLE_OPENMP "Enable OpenMP support" OFF)
option(UST_ENABLE_MPI "Enable MPI support" OFF)
option(UST_ENABLE_CUDA "Enable CUDA support" OFF)
option(UST_BUILD_LIB_HYPRE "Build HYPRE library" OFF)
option(UST_BUILD_TOOLS "Build tools" ON)
option(UST_BUILD_TESTS "Build tests" ON)

# OpenMP support
if(UST_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found")
    else()
        message(STATUS "OpenMP not found")
    endif()
else()
    message(STATUS "OpenMP disabled")
endif()

# MPI support
if(UST_ENABLE_MPI)
    find_package(MPI QUIET)
    if(MPI_CXX_FOUND)
        message(STATUS "MPI found")
    else()
        message(STATUS "MPI not found")
    endif()
else()
    message(STATUS "MPI disabled")
endif()

# Include installation directories and packaging helpers
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Subdirectories to core library and main programs
add_subdirectory(src)
add_subdirectory(applications)

# Subdirectory to additional tools
if(UST_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

# Subdirectory to unit tests
if(UST_BUILD_TESTS)
    add_subdirectory(tests)
endif()
