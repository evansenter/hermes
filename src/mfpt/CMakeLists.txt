find_package(LAPACK REQUIRED)

include_directories(c)
include_directories(h)

add_library(mfpt_static STATIC c/mfpt_params.c c/mfpt_parser.c c/mfpt_initializers.c c/mfpt_functions.c)
add_library(mfpt_shared SHARED c/mfpt_params.c c/mfpt_parser.c c/mfpt_initializers.c c/mfpt_functions.c)
add_executable(RNAmfpt c/mfpt.c c/mfpt_params.c c/mfpt_parser.c c/mfpt_initializers.c c/mfpt_functions.c)

target_link_libraries(RNAmfpt ${LAPACK_LIBRARIES})

target_link_libraries(mfpt_shared ${LAPACK_LIBRARIES})

install(TARGETS RNAmfpt DESTINATION ${HERMES_SOURCE_DIR}/bin)
install(TARGETS mfpt_static DESTINATION ${HERMES_SOURCE_DIR}/lib)
install(TARGETS mfpt_shared DESTINATION ${HERMES_SOURCE_DIR}/lib)
