find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(c)
include_directories(h)
include_directories(${FFTW_INCLUDES})
include_directories(${GSL_INCLUDES})

link_directories(${HERMES_SOURCE_DIR}/src/tpl)
link_directories(${HERMES_SOURCE_DIR}/src/mfpt)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(RNAeq c/population.c c/population_params.c c/population_initializers.c c/population_functions.c)
add_library(population STATIC c/population_params.c c/population_initializers.c c/population_functions.c)
target_link_libraries(RNAeq ${LAPACK_LIBRARIES})
target_link_libraries(RNAeq ${GSL_LIBRARIES})
target_link_libraries(RNAeq ${RNA})
target_link_libraries(RNAeq tpl)
target_link_libraries(RNAeq mfpt)

install(TARGETS RNAeq DESTINATION ${HERMES_SOURCE_DIR}/bin)
install(TARGETS population DESTINATION ${HERMES_SOURCE_DIR}/lib)
