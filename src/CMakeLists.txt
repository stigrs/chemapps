include_directories (${ChemApps_SOURCE_DIR})
include_directories (${ARMADILLO_INCLUDE_DIRS})
link_directories (${ChemApps_BINARY_DIR}/chem)

set(PROGRAMS molrot zmatxyz)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries (${program} chem ${ARMADILLO_LIBRARIES}) 
endforeach()

