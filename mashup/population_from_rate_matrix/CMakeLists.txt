find_package(OpenMP REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)

include_directories(c)
include_directories(h)
include_directories(${GSL_INCLUDES})

link_directories(${HERMES_SOURCE_DIR}/src/tpl)
link_directories(${HERMES_SOURCE_DIR}/src/multi_param)
link_directories(${HERMES_SOURCE_DIR}/src/mfpt)
link_directories(${HERMES_SOURCE_DIR}/src/population)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-write-strings")

add_executable(RateEq population_from_rate_matrix.c)
target_link_libraries(RateEq ${LAPACK_LIBRARIES})
target_link_libraries(RateEq ${GSL_LIBRARIES})
target_link_libraries(RateEq ${RNA})
target_link_libraries(RateEq tpl)
target_link_libraries(RateEq multi_param)
target_link_libraries(RateEq mfpt)
target_link_libraries(RateEq population)

install(TARGETS RateEq DESTINATION ${HERMES_SOURCE_DIR}/bin)
