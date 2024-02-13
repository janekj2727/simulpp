# include header directories
include_directories(${simul++_SOURCE_DIR}/utils/incl)
include_directories(${simul++_SOURCE_DIR}/mymath/incl)
include_directories(${simul++_SOURCE_DIR}/system/incl)

# add libraries (for each class)
add_library(utils general_utils.cpp file_openclose.cpp)
add_library(FieldFile FieldFile.cpp)
add_library(PlbFile PlbFile.cpp)
target_link_libraries(PlbFile mymath utils)
target_link_libraries(FieldFile System)

# add executable plbatomdist
add_executable(plbatomdist plbatomdist.cpp)
target_link_libraries(plbatomdist PlbFile)

# add executable mergeconfig
add_executable(mergeconfig merge_config.cpp)
target_link_libraries(mergeconfig System FieldFile)
