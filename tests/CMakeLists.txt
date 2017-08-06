include_directories (${ChemApps_SOURCE_DIR})
include_directories (${ARMADILLO_INCLUDE_DIRS})
link_directories (${ChemApps_BINARY_DIR}/chem)

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
)
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")

# Copy *.inp files to build directory.
foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach(inp_file_name)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries (${program} chem ${ARMADILLO_LIBRARIES}) 
endforeach(program)

foreach(program ${PROGRAMS})
    add_test(NAME ${program} COMMAND ${program})
endforeach(program)
