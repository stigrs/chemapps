cmake_minimum_required(VERSION 3.4.3)
project(ChemLib CXX)

include_directories(${ChemApps_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${MKL_INCLUDE_DIR})
include_directories(${SRS_INCLUDE_DIR})
link_directories(${ChemApps_BINARY_DIR})
link_directories(${SRS_LIB_DIR})

set(
    SRC_FILES
    collision.cpp
    gauss_data.cpp
    mcmm.cpp
    molecule.cpp
    molecule_io.cpp
    molrot.cpp
    molvib.cpp
    mopac.cpp
    ptable.cpp
    thermochem.cpp
    thermodata.cpp
    torsion.cpp
    troe.cpp
    tst.cpp
    tunnel.cpp
    units.cpp
    whitten_rabino.cpp
    zmatrix.cpp
)
add_library(
    chem 
    ${SRC_FILES}
)
if(MSVC)
    target_link_libraries(
        chem 
        srs.lib
    )
else()
    target_link_libraries(
        chem
        srs
    ) 
endif()

install(
    TARGETS chem
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
