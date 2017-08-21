include_directories (${ChemApps_SOURCE_DIR})
include_directories (${ARMADILLO_INCLUDE_DIRS})
link_directories (${ChemApps_BINARY_DIR})

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
add_library (
    chem 
    ${SRC_FILES}
) 
target_link_libraries (chem ${ARMADILLO_LIBRARIES}) 

install(
    TARGETS chem
    LIBRARY DESTINATION lib
)
foreach(file ${HEADER_FILES})
    install(FILES ${file} DESTINATION include/chem/)
endforeach()

