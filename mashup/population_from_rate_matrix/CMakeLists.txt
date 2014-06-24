find_package(OpenMP REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)

include_directories(c)
include_directories(h)
include_directories(${GSL_INCLUDES})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-unused-parameter -Wno-write-strings")

add_executable(RateEq population_from_rate_matrix.c)

target_link_libraries(RateEq population_static)
target_link_libraries(RateEq mfpt_static)
target_link_libraries(RateEq multi_param_static)
target_link_libraries(RateEq tpl_static)
target_link_libraries(RateEq ${LAPACK_LIBRARIES})
target_link_libraries(RateEq ${GSL_LIBRARIES})
target_link_libraries(RateEq ${RNA})
target_link_libraries(RateEq m)

install(TARGETS RateEq DESTINATION bin)
