# TODO: Install energy files from misc to bin.

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

cmake_minimum_required(VERSION 2.6)

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()

project(HERMES)

if (CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${HERMES_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${HERMES_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${HERMES_SOURCE_DIR}/lib)

include_directories(${HERMES_SOURCE_DIR}/h)

macro(c_compiler_flag_test flag_name flag fallback_flag required)
  check_c_compiler_flag(${flag} ${flag_name})
  if("${flag_name}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  elseif(NOT "${fallback_flag}" STREQUAL "")
    c_compiler_flag_test(${flag_name}_FALLBACK ${fallback_flag} "" ${required})
  elseif("${required}")
    message(FATAL_ERROR "${flag} support is required by ${CMAKE_C_COMPILER}")
  endif()
endmacro(c_compiler_flag_test)

macro(cxx_compiler_flag_test flag_name flag fallback_flag required)
  check_cxx_compiler_flag(${flag} ${flag_name})
  if("${flag_name}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  elseif(NOT "${fallback_flag}" STREQUAL "")
    cxx_compiler_flag_test(${flag_name}_FALLBACK ${fallback_flag} "" ${required})
  elseif("${required}")
    message(FATAL_ERROR "${flag} support is required by ${CMAKE_CXX_COMPILER}")
  endif()
endmacro(cxx_compiler_flag_test)

c_compiler_flag_test(C_STD_GNU99 "-std=gnu99" "" 1)
c_compiler_flag_test(C_PEDANTIC "-pedantic" "" 0)
c_compiler_flag_test(C_WALL "-Wall" "" 0)
c_compiler_flag_test(C_WEXTRA "-Wextra" "" 0)
c_compiler_flag_test(C_MARCH "-march=native" "" 0)
c_compiler_flag_test(C_CLANG_ASSEMBLER "-Wa,-q" "" 0)

cxx_compiler_flag_test(CXX_STD_CXX98 "-std=c++98" "" 1)
cxx_compiler_flag_test(CXX_PEDANTIC "-pedantic" "" 0)
cxx_compiler_flag_test(CXX_WALL "-Wall" "" 0)
cxx_compiler_flag_test(CXX_WEXTRA "-Wextra" "" 0)
cxx_compiler_flag_test(CXX_MARCH "-march=native" "" 0)
cxx_compiler_flag_test(CXX_CLANG_ASSEMBLER "-Wa,-q" "" 0)

if (CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
  c_compiler_flag_test(C_DEBUG "-g" "" 1)

  cxx_compiler_flag_test(CXX_DEBUG "-g" "" 1)
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "SAFE")
  c_compiler_flag_test(C_O2 "-O2" "" 0)

  cxx_compiler_flag_test(CXX_O2 "-O2" "" 0)
else()
  c_compiler_flag_test(C_FUNROLL_LOOPS "-funroll-loops" "" 0)
  c_compiler_flag_test(C_OFAST "-Ofast" "-O3" 1)
  c_compiler_flag_test(C_FAST_MATH "-ffast-math" "" 0)

  cxx_compiler_flag_test(CXX_FUNROLL_LOOPS "-funroll-loops" "" 0)
  cxx_compiler_flag_test(CXX_OFAST "-Ofast" "-O3" 1)
  cxx_compiler_flag_test(C_FAST_MATH "-ffast-math" "" 0)
endif()

add_subdirectory(src/tpl)
add_subdirectory(src/multi_param)
add_subdirectory(src/klp_matrix)
add_subdirectory(src/fftbor2d)
add_subdirectory(src/mfpt)
add_subdirectory(src/population)

add_subdirectory(ext/mfpt_from_fftbor2d)
add_subdirectory(ext/population_from_fftbor2d)

INSTALL(FILES misc/rna_turner1999.par misc/rna_turner2004.par DESTINATION bin)