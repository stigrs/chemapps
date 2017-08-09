include_directories (${ChemApps_SOURCE_DIR})
include_directories (${ARMADILLO_INCLUDE_DIRS})
link_directories (${ChemApps_BINARY_DIR})

set(
    SRC_FILES
    argparse.cpp 
    input.cpp
    ptable.cpp
    datum.cpp
    molecule_io.cpp
    molecule.cpp
    molrot.cpp
    math.cpp
    zmatrix.cpp
    molvib.cpp
    torsion.cpp
    mopac.cpp
)
set(
    HEADER_FILES
    argparse.h 
    input.h
    ptable.h
    datum.h
    molecule_io.h
    molecule.h
    molrot.h
    math.h
    zmatrix.h
    molvib.h
    torsion.h
    mopac.h
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

