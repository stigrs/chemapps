cmake_minimum_required(VERSION 3.4.3)
project(libchem CXX)

include_directories(${ChemApps_SOURCE_DIR}/include)
include_directories(${Stdutils_INCLUDE_DIRS})
include_directories(${Numlib_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${ChemApps_BINARY_DIR})
link_directories(${Numlib_LIBRARY_DIRS})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()

set(
    SRC_FILES
    elec_state.cpp
    geometry.cpp
    molecule.cpp
    molecule_io.cpp
    periodic_table.cpp
    rotation.cpp
    torsion.cpp
    vibration.cpp
    zmatrix.cpp
)
add_library(
    chem 
    ${SRC_FILES}
)
add_dependencies(chem Stdutils Numlib)
target_link_libraries(
    chem 
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${Numlib_LIBRARIES} 
) 

install(
    TARGETS chem
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
