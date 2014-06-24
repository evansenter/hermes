add_library(multi_param_static STATIC multi_param.c)

install(TARGETS multi_param_static DESTINATION lib)
