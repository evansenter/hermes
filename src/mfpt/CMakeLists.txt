find_package(LAPACK REQUIRED)

include_directories(c)
include_directories(h)

add_executable(RNAmfpt c/mfpt.c c/mfpt_params.c c/mfpt_parser.c c/mfpt_initializers.c c/mfpt_functions.c)
add_library(mfpt STATIC c/mfpt_params.c c/mfpt_parser.c c/mfpt_initializers.c c/mfpt_functions.c)
target_link_libraries(RNAmfpt ${LAPACK_LIBRARIES})

install(TARGETS RNAmfpt DESTINATION ${HERMES_SOURCE_DIR}/bin)
install(TARGETS mfpt DESTINATION ${HERMES_SOURCE_DIR}/lib)
