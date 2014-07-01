find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_package(GSL REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(${FFTW_INCLUDES})
include_directories(${GSL_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-write-strings")

add_executable(FFTmfpt mfpt_from_fftbor2d.cpp)

target_link_libraries(FFTmfpt fftbor_static)
target_link_libraries(FFTmfpt mfpt_static)
target_link_libraries(FFTmfpt multi_param_static)
target_link_libraries(FFTmfpt ${FFTW_LIBRARIES})
target_link_libraries(FFTmfpt ${GSL_LIBRARIES})
target_link_libraries(FFTmfpt ${RNA})
target_link_libraries(FFTmfpt m)

install(TARGETS FFTmfpt DESTINATION bin)
