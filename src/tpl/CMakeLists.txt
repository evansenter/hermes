set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare")

add_library(tpl_static STATIC tpl.c)

install(TARGETS tpl_static DESTINATION lib)
