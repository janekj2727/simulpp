# include header directories
include_directories(
    ${simul++_SOURCE_DIR}/system/incl
    ${simul++_SOURCE_DIR}/mymath/incl
    ${simul++_SOURCE_DIR}/utils/incl
)

# add library System
add_library(System Atom.cpp Molecule.cpp HarmonicBond.cpp HarmonicAngle.cpp CosineDihedral.cpp ConstraintBond.cpp SimulatedSystem.cpp Measurement.cpp AllPairList.cpp VerletList.cpp Pair.cpp LJsystem.cpp CutoffElstat.cpp EwaldElstat.cpp LJBond.cpp ElstatBond.cpp EwaldSplinesElstat.cpp MeanSquareDispl.cpp)
target_link_libraries(System mymath utils FieldFile)

# # add libraries (for each class)
# add_library(Atom Atom.cpp)
# target_link_libraries(Atom Vector Matrix)
# add_library(Bond Bond.cpp)
# target_link_libraries(Bond Molecule)
# add_library(ConstraintBond ConstraintBond.cpp)
# target_link_libraries(ConstraintBond Molecule)
# add_library(LJsystem LJsystem.cpp)
# target_link_libraries(LJsystem NewtonMethod)
# add_library(Molecule Molecule.cpp)
# target_link_libraries(Molecule Atom)
# add_library(SimulatedSystem SimulatedSystem.cpp)
# target_link_libraries(SimulatedSystem Molecule LJsystem math_utils file_openclose)