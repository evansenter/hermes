find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(cpp)
include_directories(h)
include_directories(${FFTW_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(fftbor_static STATIC cpp/fftbor2d_params.cpp cpp/fftbor2d_initializers.cpp cpp/fftbor2d_functions.cpp)
add_library(fftbor_shared SHARED cpp/fftbor2d_params.cpp cpp/fftbor2d_initializers.cpp cpp/fftbor2d_functions.cpp)
add_executable(FFTbor2D cpp/fftbor2d.cpp cpp/fftbor2d_params.cpp cpp/fftbor2d_initializers.cpp cpp/fftbor2d_functions.cpp)

target_link_libraries(FFTbor2D ${FFTW_LIBRARIES})
target_link_libraries(FFTbor2D ${RNA})

target_link_libraries(fftbor_shared ${FFTW_LIBRARIES})
target_link_libraries(fftbor_shared ${RNA})

install(TARGETS FFTbor2D DESTINATION bin)
install(TARGETS fftbor_static DESTINATION lib)
install(TARGETS fftbor_shared DESTINATION lib)
