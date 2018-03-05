cmake_minimum_required(VERSION 3.4.3)
project(ChemProgs CXX)

include_directories(${ChemApps_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIR})
include_directories(${SRS_INCLUDE_DIR})
link_directories(${ChemApps_BINARY_DIR}/chem)
link_directories(${MKL_LIB_DIR})
link_directories(${SRS_LIB_DIR})
link_directories($ENV{BOOST_LIBRARYDIR})

set(
    PROGRAMS 
    brw
    mcmm 
    molrot 
    thermo
    tst
    zmatxyz
)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    if (WIN32)
        target_link_libraries (
            ${program} 
            chem 
            srs.lib
            mkl_intel_c.lib 
            mkl_sequential.lib 
            mkl_core.lib
            ${Boost_LIBRARIES}
        ) 
    else()
        target_link_libraries (
            ${program} 
            chem 
            srs
            ${Boost_LIBRARIES}
        )
    endif()
endforeach()

install(
    TARGETS ${PROGRAMS}
    RUNTIME DESTINATION bin
    CONFIGURATIONS Release
)
