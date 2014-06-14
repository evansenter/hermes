SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare")

add_library(tpl STATIC tpl.c)
