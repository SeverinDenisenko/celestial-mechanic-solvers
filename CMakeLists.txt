cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(AdamsCelestialMechanics VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost REQUIRED)
find_package(MPFR REQUIRED)

set(ADAMS_SOLVER_LIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/lib/adams.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lib/rk.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/lib/adams_mpfr.cpp)

add_library(AdamsSolverLib STATIC ${ADAMS_SOLVER_LIB_SRC})
target_include_directories(AdamsSolverLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/ ${Boost_INCLUDE_DIR})
target_link_libraries(AdamsSolverLib ${MPFR_LIBRARIES} ${Boost_LIBRARIES})
 target_compile_options(AdamsSolverLib PRIVATE -Wall -Wextra -Wpedantic -Werror -Wno-unused-private-field -Wno-deprecated-declarations -Wno-deprecated-copy-with-user-provided-copy)

set(ADAMS_SOLVER_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
add_executable(AdamsSolver ${ADAMS_SOLVER_SRC})

target_link_libraries(AdamsSolver AdamsSolverLib)
