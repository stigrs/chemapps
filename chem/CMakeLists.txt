cmake_minimum_required(VERSION 3.4.3)
project(ChemLib CXX)

include_directories(${ChemApps_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${SRS_INCLUDE_DIRS})
link_directories(${ChemApps_BINARY_DIR})

set(
    SRC_FILES
    #collision.cpp
    ptable.cpp
    #mcmm.cpp
    #molecule.cpp
    molecule_io.cpp
    #molrot.cpp
    #molvib.cpp
    #mopac.cpp
    #thermochem.cpp
    #thermodata.cpp
    #torsion.cpp
    #tst.cpp
    #tunnel.cpp
    #zmatrix.cpp
)
add_library(
    chem 
    ${SRC_FILES}
)
if(MSVC)
    target_link_libraries(
        chem 
    )
else()
    target_link_libraries (chem) 
endif()

install(
    TARGETS chem
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    CONFIGURATIONS Release
)
