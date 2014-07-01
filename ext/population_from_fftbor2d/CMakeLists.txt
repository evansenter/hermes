find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_package(GSL REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(${FFTW_INCLUDES})
include_directories(${GSL_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-write-strings")

add_executable(FFTeq population_from_fftbor2d.cpp)

target_link_libraries(FFTeq fftbor_static)
target_link_libraries(FFTeq population_static)
target_link_libraries(FFTeq mfpt_static)
target_link_libraries(FFTeq multi_param_static)
target_link_libraries(FFTeq tpl_static)
target_link_libraries(FFTeq ${FFTW_LIBRARIES})
target_link_libraries(FFTeq ${GSL_LIBRARIES})
target_link_libraries(FFTeq ${RNA})
target_link_libraries(FFTeq m)

install(TARGETS FFTeq DESTINATION bin)
