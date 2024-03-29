# include header directories
include_directories(
    ${simul++_SOURCE_DIR}/integrators/incl
    ${simul++_SOURCE_DIR}/mymath/incl
    ${simul++_SOURCE_DIR}/system/incl
    ${simul++_SOURCE_DIR}/utils/incl
)

# add libraries (for each class) and specify their dependencies
add_library(AbstractGearIntegrator AbstractGearIntegrator.cpp)
target_link_libraries(AbstractGearIntegrator mymath System)
add_library(AbstractVerletIntegrator AbstractVerletIntegrator.cpp)
target_link_libraries(AbstractVerletIntegrator mymath System)
add_library(GearIntegrator GearIntegrator.cpp)
target_link_libraries(GearIntegrator AbstractGearIntegrator)
add_library(GearNTVinitBerendsen GearNTVinitBerendsen.cpp)
target_link_libraries(GearNTVinitBerendsen AbstractGearIntegrator)
add_library(GearNVTBerendsen GearNVTBerendsen.cpp)
target_link_libraries(GearNVTBerendsen AbstractGearIntegrator)
add_library(GearNVTNose GearNVTNose.cpp)
target_link_libraries(GearNVTNose AbstractGearIntegrator)
add_library(GearNPTBerendsen GearNPTBerendsen.cpp)
target_link_libraries(GearNPTBerendsen AbstractGearIntegrator)
add_library(GearNPTNose GearNPTNose.cpp)
target_link_libraries(GearNPTNose AbstractGearIntegrator)
add_library(VerletIntegrator VerletIntegrator.cpp)
target_link_libraries(VerletIntegrator AbstractVerletIntegrator)
add_library(VerletNTVinitBerendsen VerletNTVinitBerendsen.cpp)
target_link_libraries(VerletNTVinitBerendsen AbstractVerletIntegrator)
add_library(VerletNVTBerendsen VerletNVTBerendsen.cpp)
target_link_libraries(VerletNVTBerendsen AbstractVerletIntegrator)
add_library(VerletNVTNose VerletNVTNose.cpp)
target_link_libraries(VerletNVTNose AbstractVerletIntegrator)
add_library(VerletNVTNoseTRVP VerletNVTNoseTRVP.cpp)
target_link_libraries(VerletNVTNoseTRVP AbstractVerletIntegrator)
add_library(VerletNPTBerendsen VerletNPTBerendsen.cpp)
target_link_libraries(VerletNPTBerendsen AbstractVerletIntegrator)
add_library(VerletNPTNose VerletNPTNose.cpp)
target_link_libraries(VerletNPTNose AbstractVerletIntegrator)
add_library(VerletNPTNoseTRVP VerletNPTNoseTRVP.cpp)
target_link_libraries(VerletNPTNoseTRVP AbstractVerletIntegrator)
add_library(VerletNVTNoseIter VerletNVTNoseIter.cpp)
target_link_libraries(VerletNVTNoseIter AbstractVerletIntegrator)
add_library(VerletNPTNoseIter VerletNPTNoseIter.cpp)
target_link_libraries(VerletNPTNoseIter AbstractVerletIntegrator)
add_library(VerletNPTAlejandre VerletNPTAlejandre.cpp)
target_link_libraries(VerletNPTAlejandre AbstractVerletIntegrator)
