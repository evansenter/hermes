find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_package(LAPACK REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(${FFTW_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-write-strings")

add_executable(FFTmfpt mfpt_from_fftbor2d.cpp)

target_link_libraries(FFTmfpt multi_param_shared)
target_link_libraries(FFTmfpt fftbor_shared)
target_link_libraries(FFTmfpt mfpt_shared)
target_link_libraries(FFTmfpt ${FFTW_LIBRARIES})
target_link_libraries(FFTmfpt ${LAPACK_LIBRARIES})
target_link_libraries(FFTmfpt ${RNA})

install(TARGETS FFTmfpt DESTINATION bin)
