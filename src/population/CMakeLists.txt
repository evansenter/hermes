find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(c)
include_directories(h)
include_directories(${FFTW_INCLUDES})
include_directories(${GSL_INCLUDES})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_library(population_static STATIC c/population_params.c c/population_initializers.c c/population_functions.c)
add_executable(RNAeq c/population.c c/population_params.c c/population_initializers.c c/population_functions.c)

target_link_libraries(RNAeq mfpt_static)
target_link_libraries(RNAeq tpl_static)
target_link_libraries(RNAeq ${LAPACK_LIBRARIES})
target_link_libraries(RNAeq ${GSL_LIBRARIES})
target_link_libraries(RNAeq ${RNA})
target_link_libraries(RNAeq m)

install(TARGETS RNAeq DESTINATION bin)
install(TARGETS population_static DESTINATION lib)
