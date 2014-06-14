cmake_minimum_required(VERSION 2.6)

SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

project(HERMES)

execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GXX_VERSION)

if (GCC_VERSION VERSION_LESS 4.8)
  message(FATAL_ERROR "gcc version < 4.8")
endif()

if (GXX_VERSION VERSION_LESS 4.8)
  message(FATAL_ERROR "g++ version < 4.8")
endif()

enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -pedantic -funroll-loops -Wall -Wextra -Ofast -march=native -Wa,-q")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -funroll-loops -Wall -Wextra -Ofast -march=native -Wa,-q")

add_subdirectory(src/tpl)
add_subdirectory(src/multi_param)
add_subdirectory(src/fftbor2d)
add_subdirectory(src/mfpt)
add_subdirectory(src/population)

add_subdirectory(mashup/mfpt_from_fftbor2d)
add_subdirectory(mashup/population_from_fftbor2d)
add_subdirectory(mashup/population_from_rate_matrix)
