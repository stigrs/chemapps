include_directories (${ChemApps_SOURCE_DIR})
include_directories (${ARMADILLO_INCLUDE_DIRS})
link_directories (${ChemApps_BINARY_DIR}/chem)

# Copy *.inp files to build directory.
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")

foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

# Prepare "Catch" library for other executables.
set(CATCH_INCLUDE_DIR /usr/local/include/catch)
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
    target_link_libraries(
        ${name}
        chem
        ${ARMADILLO_LIBRARIES}
        test_catch
    )
    add_dependencies(${name} Catch)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
    test_datum
    test_input
    test_molecule
    test_molecule_io
    test_molrot
    test_molvib
    test_pdist_matrix
    test_ptable
    test_zmatrix
    test_torsion
    test_mopac
    test_mcmm
)

foreach(program ${PROGRAMS})
    add_chemapps_test(${program})
endforeach()

