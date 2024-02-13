# include header directories
include_directories(${simul++_SOURCE_DIR}/mymath/incl)
include_directories(${simul++_SOURCE_DIR}/utils/incl)

# add libraries (for each class) and their dependencies
add_library(mymath Matrix.cpp Vector.cpp LinearSystem.cpp NewtonMethod.cpp math_utils.cpp)
add_library(ARMA ARMA.cpp)
target_link_libraries(ARMA mymath)

add_library(interp LinearInterpolator.cpp Linear3PInterpolator.cpp Linear4PInterpolator.cpp HermiteCubSplines.cpp NaturalCubSplines.cpp MacsimusQuadSplines.cpp MacsimusHyperbSplines.cpp)

# add executables
add_executable(arma arma_main.cpp)
target_link_libraries(arma ARMA)
target_link_libraries(arma utils)

add_executable(newton test_newton.cpp)
target_link_libraries(newton mymath)

add_executable(mymath_script EXCLUDE_FROM_ALL mymath_main.cpp)
target_link_libraries(mymath_script ARMA interp)

# add_executable(newton test_newton.cpp)
# target_link_libraries(newton mymath)

# add_executable(test_rand test_rand.cpp)
# target_link_libraries(test_rand mymath)

# # add libraries (for each class) and their dependencies
# add_library(ARMA ARMA.cpp)
# target_link_libraries(ARMA Matrix LinearSystem Vector math_utils)
# add_library(LinearSystem LinearSystem.cpp)
# target_link_libraries(LinearSystem Vector Matrix)
# add_library(math_utils math_utils.cpp)
# add_library(Matrix Matrix.cpp)
# target_link_libraries(Matrix Vector)
# add_library(NewtonMethod NewtonMethod.cpp)
# target_link_libraries(NewtonMethod Matrix Vector LinearSystem)
# add_library(Vector Vector.cpp)