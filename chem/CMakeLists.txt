cmake_minimum_required(VERSION 3.4.3)
project(ChemLib CXX)

include_directories(${ChemApps_SOURCE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIR})
link_directories(${ChemApps_BINARY_DIR})
if(WIN32)
    link_directories($ENV{MKLROOT})
endif()

set(
    SRC_FILES
    datum.cpp
    input.cpp
    ptable.cpp
    math.cpp
    mcmm.cpp
    molecule.cpp
    molecule_io.cpp
    molrot.cpp
    molvib.cpp
    mopac.cpp
    thermochem.cpp
    thermodata.cpp
    torsion.cpp
    zmatrix.cpp
)
set(
    HEADER_FILES
    arma_io.h
    conformer.h
    datum.h
    input.h
    math.h
    mcmm.h
    molecule.h
    molecule_io.h
    molrot.h
    molvib.h
    mopac.h
    ptable.h
    thermodata.h
    torsion.h
    zmatrix.h
)
add_library(
    chem 
    ${SRC_FILES}
) 
if(WIN32)
    target_link_libraries(
        chem 
        mkl_intel_c.lib
        mkl_sequential.lib
        mkl_core.lib
    )
else()
    target_link_libraries (chem ${ARMADILLO_LIBRARIES}) 
endif()

install(
    TARGETS chem
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
foreach(file ${HEADER_FILES})
    install(FILES ${file} DESTINATION include/chem/)
endforeach()

