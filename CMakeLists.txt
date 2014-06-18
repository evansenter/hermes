include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

cmake_minimum_required(VERSION 2.6)

SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

project(HERMES)

execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GXX_VERSION)

if(GCC_VERSION VERSION_LESS 4.8)
  message(FATAL_ERROR "gcc version < 4.8")
endif()

if(GXX_VERSION VERSION_LESS 4.8)
  message(FATAL_ERROR "g++ version < 4.8")
endif()

enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

foreach(include_path ${ADDITIONAL_INCLUDE_PATHS})
  message(STATUS "Adding include path: ${include_path}")
  include_directories(${include_path})
endforeach(include_path ${ADDITIONAL_INCLUDE_PATHS})

foreach(library_path ${ADDITIONAL_LIBRARY_PATHS})
  message(STATUS "Adding library path: ${library_path}")
  link_directories(${library_path})
endforeach(library_path ${ADDITIONAL_LIBRARY_PATHS})

macro(c_compiler_flag_test flag_name flag fallback_flag)
  check_c_compiler_flag(${flag} ${flag_name})
  if(${flag_name})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  elseif(${fallback_flag})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${fallback_flag}")
  endif()
endmacro(c_compiler_flag_test)

macro(cxx_compiler_flag_test flag_name flag fallback_flag)
  check_cxx_compiler_flag(${flag} ${flag_name})
  if(${flag_name})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  elseif(${fallback_flag})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${fallback_flag}")
  endif()
endmacro(cxx_compiler_flag_test)

c_compiler_flag_test(C_STD_GNU99 "-std=gnu99" "")
c_compiler_flag_test(C_PEDANTIC "-pedantic" "")
c_compiler_flag_test(C_FUNROLL_LOOPS "-funroll-loops" "")
c_compiler_flag_test(C_WALL "-Wall" "")
c_compiler_flag_test(C_WEXTRA "-Wextra" "")
c_compiler_flag_test(C_OFAST "-Ofast" "-O3")
c_compiler_flag_test(C_MARCH "-march=native" "")
c_compiler_flag_test(C_CLANG_ASSEMBLER "-Wa,-q" "")

cxx_compiler_flag_test(CXX_STD_CXX11 "-std=c++11" "")
cxx_compiler_flag_test(CXX_PEDANTIC "-pedantic" "")
cxx_compiler_flag_test(CXX_FUNROLL_LOOPS "-funroll-loops" "")
cxx_compiler_flag_test(CXX_WALL "-Wall" "")
cxx_compiler_flag_test(CXX_WEXTRA "-Wextra" "")
cxx_compiler_flag_test(CXX_OFAST "-Ofast" "-O3")
cxx_compiler_flag_test(CXX_MARCH "-march=native" "")
cxx_compiler_flag_test(CXX_CLANG_ASSEMBLER "-Wa,-q" "")

add_subdirectory(src/tpl)
add_subdirectory(src/multi_param)
add_subdirectory(src/fftbor2d)
add_subdirectory(src/mfpt)
add_subdirectory(src/population)

add_subdirectory(mashup/mfpt_from_fftbor2d)
add_subdirectory(mashup/population_from_fftbor2d)
add_subdirectory(mashup/population_from_rate_matrix)
