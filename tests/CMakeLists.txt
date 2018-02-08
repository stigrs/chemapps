cmake_minimum_required(VERSION 3.4.3)
project(ChemTests CXX)

include_directories (${ChemApps_SOURCE_DIR})
include_directories (${MKL_INCLUDE_DIR})
include_directories (${SRS_INCLUDE_DIR})
link_directories (${ChemApps_BINARY_DIR}/chem)
link_directories (${MKL_LIB_DIR})
link_directories (${SRS_LIB_DIR})

# Copy *.inp files to build directory.
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")

foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

# Prepare "Catch" library for other executables.
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Create the main executable for each test.
add_library(test_catch STATIC test.cpp)
target_link_libraries(
    test_catch
    Catch
)

function(add_chemapps_test name)
    add_executable(${name} ${name}.cpp)
    if (MSVC)
        target_link_libraries (
            ${program} 
            chem 
            srs.lib
            test_catch
        ) 
    else()
        target_link_libraries (
            ${program} 
            chem 
            srs
            test_catch
        )
    endif()
    add_dependencies(${name} Catch)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
    #test_mcmm
    test_molecule
    test_molecule_io
    test_molrot
    test_molvib
    #test_mopac
    test_ptable
    #test_thermochem
    #test_thermodata
    test_torsion
    #test_tst
    #test_tunnel
    test_zmatrix
)

foreach(program ${PROGRAMS})
    add_chemapps_test(${program})
endforeach()
