find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)
find_package(LAPACK REQUIRED)
find_library(RNA libRNA.a PATHS)

include_directories(${FFTW_INCLUDES})

# link_directories(${HERMES_SOURCE_DIR}/src/multi_param)
# link_directories(${HERMES_SOURCE_DIR}/src/fftbor2d)
# link_directories(${HERMES_SOURCE_DIR}/src/mfpt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-write-strings")

add_executable(FFTmfpt mfpt_from_fftbor2d.cpp)
target_link_libraries(FFTmfpt multi_param)
target_link_libraries(FFTmfpt fftbor)
target_link_libraries(FFTmfpt mfpt)
target_link_libraries(FFTmfpt ${FFTW_LIBRARIES})
target_link_libraries(FFTmfpt ${LAPACK_LIBRARIES})
target_link_libraries(FFTmfpt ${RNA})

install(TARGETS FFTmfpt DESTINATION ${HERMES_SOURCE_DIR}/bin)
