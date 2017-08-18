cmake_minimum_required(VERSION 3.4.3)
project(ChemApps CXX)

# By default, ChemApps is install in the HOME directory.
set(CMAKE_INSTALL_PREFIX $ENV{HOME})

# Export compile commands.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enforce C++1 standard.
set(CMAKE_CXX_STANDARD 11)

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
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pedantic -Wall -Wshadow")
        add_definitions(-DARMA_NO_DEBUG)
    endif()
    if(APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -framework Accelerate")
    endif()
endif()

# Set debugging compiler options.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -pedantic -Wall -Wshadow")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
        add_definitions(-DARMA_EXTRA_DEBUG)
    endif()
endif()

# Armadillo library is required.
find_package(Armadillo 7.950.1 REQUIRED)

# Boost library is required.
find_package(Boost 1.58.0 COMPONENTS program_options REQUIRED)
set(Boost_USE_STATIC_LIBS OFF)

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

