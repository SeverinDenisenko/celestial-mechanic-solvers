cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(AdamsCelestialMechanics VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost REQUIRED)
find_package(MPFR REQUIRED)

set(ODE_SOLVER_LIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/lib/adams_extrapolation_solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/adams_interpolation_solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/adams_predictor_corrector_solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/jacoby_matrix_evaluators.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/full_pivot_gauss_solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/newton_root_finder.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lib/rk4_solver.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lib/quad_integrator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lib/simpson_integrator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lib/utils.cpp)

add_compile_definitions("-DSNUMLIB_USE_MPFR=1")

add_library(ODESolversLib STATIC ${ODE_SOLVER_LIB_SRC})
target_include_directories(ODESolversLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/ ${CMAKE_CURRENT_SOURCE_DIR}/deps/ ${Boost_INCLUDE_DIR})
target_link_libraries(ODESolversLib PUBLIC ${MPFR_LIBRARIES} ${Boost_LIBRARIES})
target_compile_options(ODESolversLib PUBLIC -Wall -Wextra -Wpedantic -Wfloat-conversion -Wno-unused-private-field -Wno-deprecated-declarations -Wno-deprecated-copy-with-user-provided-copy -O3)
target_link_options(ODESolversLib PUBLIC)

set(SOLVER_2_BODY_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/solver_2_body.cpp)
add_executable(Solver2Body ${SOLVER_2_BODY_SRC})
target_compile_options(Solver2Body PRIVATE -Wall -Wextra -O3)
target_link_libraries(Solver2Body ODESolversLib)

add_executable(matrix_solver_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/matrix_solver_test.cpp)
target_compile_options(matrix_solver_test PRIVATE -Wall -Wextra -O3)
target_link_libraries(matrix_solver_test ODESolversLib)

add_executable(root_finder_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/root_finder_test.cpp)
target_compile_options(root_finder_test PRIVATE -Wall -Wextra -O3)
target_link_libraries(root_finder_test ODESolversLib)

set(SOLVER_N_BODY_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/solver_n_body.cpp)
add_executable(SolverNBody ${SOLVER_N_BODY_SRC})
target_compile_options(SolverNBody PRIVATE -Wall -Wextra  -O3)
target_link_libraries(SolverNBody ODESolversLib)
