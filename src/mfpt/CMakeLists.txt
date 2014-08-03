find_package(GSL REQUIRED)

include_directories(c)
include_directories(h)
include_directories(${GSL_INCLUDES})

add_library(mfpt_static STATIC c/mfpt_params.c c/mfpt_functions.c)
add_executable(RNAmfpt c/mfpt.c c/mfpt_params.c c/mfpt_functions.c)

target_link_libraries(RNAmfpt klp_matrix_static)
target_link_libraries(RNAmfpt ${GSL_LIBRARIES})
target_link_libraries(RNAmfpt m)

install(TARGETS RNAmfpt DESTINATION bin)
install(TARGETS mfpt_static DESTINATION lib)
