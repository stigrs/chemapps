cmake_minimum_required(VERSION 3.4.3)
project(ChemApps CXX)

# By default, ChemApps is install in the HOME directory.
set(CMAKE_INSTALL_PREFIX $ENV{HOME})

# Export compile commands.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enforce C++14 standard (needed for GSL).
set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(MSVC)
    option(BUILD_SHARED_LIBS "Compile shared libraries." OFF)
else()
    option(BUILD_SHARED_LIBS "Compile shared libraries." ON)
endif()
option(BUILD_TESTS "Build tests." ON)

# Set release compiler options.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic -Wall -Wshadow -Wextra -Werror")
    endif()
    if(APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -framework Accelerate")
    endif()
    if(MSVC)
        set(CMAKE_CXX_FLAGS_RELEASE "/std:c++14 /W4 /WX /GR /EHsc /MT /O2")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/force:multiple /ignore:4075")
    endif()
    add_definitions(-DNDEBUG)
    add_definitions(-DARMA_NO_DEBUG)
    add_definitions(-DGSL_THROW_ON_CONTRACT_VIOLATION)
endif()

# Set debugging compiler options.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -Wpedantic -Wall -Wshadow -Wextra -Werror")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
    endif()
    if(MSVC)
        set(CMAKE_CXX_FLAGS_DEBUG "/std:c++14 /EHsc /W4 /WX /MTd /GR /Zi /Od /RTC1")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/force:multiple /ignore:4075")
    endif()
    add_definitions(-DGSL_THROW_ON_CONTRACT_VIOLATION)
endif()

# Armadillo library is required.
if(MSVC)
    find_path(ARMADILLO_INCLUDE_DIRS armadillo HINTS $ENV{HOME}/include)
    add_definitions(-DARMA_USE_BLAS)
    add_definitions(-DARMA_USE_LAPACK)
else()
    find_package(Armadillo 7.950.1 REQUIRED)
endif()

# Boost library is required.
if(MSVC)
    set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.58.0 COMPONENTS program_options REQUIRED)

# GSL library is required.
find_path(GSL_INCLUDE_DIR gsl HINTS $ENV{HOME}/include /usr/include /usr/local/include)

add_subdirectory(chem)
add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Print status messages.
message(STATUS "CMAKE_BUILD_TYPE  = ${CMAKE_BUILD_TYPE}")
message(STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}")
message(STATUS "BUILD_TESTS       = ${BUILD_TESTS}")

