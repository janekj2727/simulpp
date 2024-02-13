# Add links to include files
include_directories (
    ${simul++_SOURCE_DIR}/main/incl
    ${simul++_SOURCE_DIR}/mymath/incl
    ${simul++_SOURCE_DIR}/integrators/incl
    ${simul++_SOURCE_DIR}/system/incl
    ${simul++_SOURCE_DIR}/utils/incl
)

# Link directories 
# link_directories(
#     ${simul++_SOURCE_DIR}/mymath/src
#     ${simul++_SOURCE_DIR}/integrators/src
#     ${simul++_SOURCE_DIR}/system/src 
#     ${simul++_SOURCE_DIR}/utils/src 
# )

# Add library Simulation
add_library(Simulation Simulation.cpp)
target_link_libraries(Simulation utils
    VerletIntegrator VerletNTVinitBerendsen VerletNVTBerendsen
    GearIntegrator GearNTVinitBerendsen GearNVTBerendsen GearNVTNose
    VerletNVTNose VerletNVTNoseTRVP VerletNPTBerendsen GearNPTBerendsen
    VerletNPTNose VerletNPTNoseTRVP GearNPTNose VerletNVTNoseIter
    VerletNPTNoseIter VerletNPTAlejandre)

# Defining executable
add_executable(simul++ simul++.cpp)

# Libraries for executable
target_link_libraries(simul++ Simulation System)

# # Libraries for executable
# target_link_libraries(simul++ Vector)
# target_link_libraries(simul++ Matrix)
# target_link_libraries(simul++ LinearSystem)
# target_link_libraries(simul++ NewtonMethod)
# target_link_libraries(simul++ math_utils)
# target_link_libraries(simul++ Atom)
# target_link_libraries(simul++ Bond)
# target_link_libraries(simul++ ConstraintBond)
# target_link_libraries(simul++ LJsystem)
# target_link_libraries(simul++ Molecule)
# target_link_libraries(simul++ SimulatedSystem)
# target_link_libraries(simul++ file_openclose)
# target_link_libraries(simul++ general_utils)
# target_link_libraries(simul++ AbstractGearIntegrator)
# target_link_libraries(simul++ AbstractVerletIntegrator)
# target_link_libraries(simul++ GearIntegrator)
# target_link_libraries(simul++ GearNTVinitBerendsen)
# target_link_libraries(simul++ GearNVTBerendsen)
# target_link_libraries(simul++ VerletIntegrator)
# target_link_libraries(simul++ VerletNTVinitBerendsen)
# target_link_libraries(simul++ VerletNVTBerendsen)
# target_link_libraries(simul++ Simulation)
