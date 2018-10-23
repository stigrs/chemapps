cmake_minimum_required(VERSION 3.4.3)
project(ChemApps CXX)

#-------------------------------------------------------------------------------
#
# Set compiler options:

# Enforce C++14 standard.
set(CMAKE_CXX_STANDARD 14)

# Export compile commands.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WIN32)
    option(BUILD_SHARED_LIBS "Compile shared libraries." OFF)
else()
    option(BUILD_SHARED_LIBS "Compile shared libraries." ON)
endif()

# Set default MSVC compiler options to avoid D9025 error.
if(MSVC)
    set(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS")
endif()

# Set release compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/DNDEBUG /std:c++14 /W4 /WX /wd4100 /GR /EHsc /MT /O2")
endif()

# Set debugging compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/std:c++14 /EHsc /W4 /WX /wd4100 /MT /GR /Zi /Od")
endif()

#-------------------------------------------------------------------------------
#
# Include external dependencies:

include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Git)

# BLAS and LAPACK is required.
if(DEFINED ENV{MKLROOT})
	find_package(MKL)
	if(MKL_FOUND)
		add_definitions(${MKL_DEFINITIONS})
		set(BLAS_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
		set(BLAS_LIBRARIES ${MKL_LIBRARIES})
		set(LAPACK_LIBRARIES "")
	endif()
else()
	set(BLA_VENDOR "OpenBLAS")
	if(WIN32)
		set(BLAS_DIR C:/local/OpenBLAS.0.2.14.1/lib/native)
		if(${CMAKE_CL_64} EQUAL 1)
			find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/x64)
			set(BLAS_BIN "${BLAS_DIR}/bin/x64")
		else()
			find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/win32)
			set(BLAS_BIN "${BLAS_DIR}/bin/win32")
		endif()
	elseif(APPLE)
		set(BLAS_DIR /usr/local/opt/openblas)
		find_file(BLAS_LIBRARIES libopenblas.dylib HINTS ${BLAS_DIR}/lib)
	else()
		# Workaround: LAPACKE is not a part of OpenBLAS on Ubuntu
		set(BLAS_DIR /usr)
		find_package(BLAS)
		find_file(LAPACK_LIBRARIES liblapacke.so HINTS ${BLAS_DIR}/lib ${BLAS_DIR}/lib/x86_64-linux-gnu) 
		find_path(BLAS_INCLUDE_DIRS cblas.h HINTS ${BLAS_DIR}/include /usr/include /usr/local/include)
		find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${BLAS_DIR}/include /usr/include /usr/local/include)
	endif()
	find_path(BLAS_INCLUDE_DIRS cblas.h HINTS ${BLAS_DIR}/include)
	find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${BLAS_DIR}/include)
endif()

# Stdutils is required.
find_package(Stdutils)
if(Stdutils_FOUND)
    add_custom_target(Stdutils)
else()
    list(APPEND STDUTILS_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/external"
        "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
        "-DBUILD_TESTS=OFF"
    )
    ExternalProject_Add(
        Stdutils
        PREFIX ${CMAKE_BINARY_DIR}/stdutils
        GIT_REPOSITORY https://github.com/stigrs/stdutils.git
        CMAKE_ARGS ${STDUTILS_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 0
    )
    set(Stdutils_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/include)
endif()

# Numlib is required.
find_package(Numlib)
if(Numlib_FOUND)
    add_custom_target(Numlib)
else()
    list(APPEND NUMLIB_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/external"
        "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
        "-DBUILD_TESTS=OFF"
        "-DBUILD_BENCH=OFF"
        "-DBUILD_EXAMPLES=OFF"
    )
    ExternalProject_Add(
        Numlib
        PREFIX ${CMAKE_BINARY_DIR}/numlib
        GIT_REPOSITORY https://github.com/stigrs/numlib.git
        CMAKE_ARGS ${NUMLIB_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 0
    )
    set(Numlib_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/include)
    set(Numlib_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/external/lib)
    set(Numlib_LIBRARIES "num")
endif()

#-------------------------------------------------------------------------------
#
# Build the software:

# Create chemlib library.
add_subdirectory(lib)

# Enable testing.
option(BUILD_TESTS "Build tests." ON)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Install target.
install(
    DIRECTORY include/chem
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

