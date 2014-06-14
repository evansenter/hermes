add_library(multi_param STATIC multi_param.c)

install(TARGETS multi_param DESTINATION ${HERMES_SOURCE_DIR}/lib)
