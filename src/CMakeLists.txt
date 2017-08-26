cmake_minimum_required(VERSION 3.4.3)
project(ChemProgs CXX)

include_directories(${ChemApps_SOURCE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${ChemApps_BINARY_DIR}/chem)
if (MSVC)
    link_directories($ENV{MKLROOT})
    link_directories($ENV{BOOST_LIBRARYDIR})
endif()

set(
    PROGRAMS 
    mcmm 
    molrot 
    zmatxyz
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    if (MSVC)
        target_link_libraries (
            ${program} 
            chem 
            mkl_intel_c.lib 
            mkl_sequential.lib 
            mkl_core.lib
            ${Boost_LIBRARIES}
        ) 
    else()
        target_link_libraries (
            ${program} 
            chem 
            ${ARMADILLO_LIBRARIES} 
            ${Boost_LIBRARIES}
        )
    endif()
endforeach()
