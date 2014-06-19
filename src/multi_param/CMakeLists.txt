add_library(multi_param_static STATIC multi_param.c)
add_library(multi_param_shared SHARED multi_param.c)

install(TARGETS multi_param_static DESTINATION lib)
install(TARGETS multi_param_shared DESTINATION lib)
