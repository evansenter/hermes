set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare")

add_library(tpl_static STATIC tpl.c)
add_library(tpl_shared SHARED tpl.c)

install(TARGETS tpl_static DESTINATION lib)
install(TARGETS tpl_shared DESTINATION lib)
