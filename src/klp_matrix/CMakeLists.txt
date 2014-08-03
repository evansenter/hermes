include_directories(c)
include_directories(h)

add_library(klp_matrix_static STATIC c/klp_matrix_parser.c c/klp_matrix_initializers.c c/klp_matrix_functions.c)

install(TARGETS klp_matrix_static DESTINATION lib)
