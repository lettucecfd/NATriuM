#
# CMake File for NATriuM test
#    

#
# It is good practise to specify a version requirement:
#
cmake_minimum_required( VERSION 2.8.8 )


# 
#
set(    TARGET 
        fft
        )


  
#
# Specify a target:
#
add_executable( ${TARGET} 
                FFT.cpp
                )
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
#
# Find FFT library:
#
find_library( FFTW_LIBRARY NAMES fftw3 fftw)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
if(UNIX AND NOT WIN32)
     find_library(FFTW_libm_LIBRARY
         NAMES m)
     list(APPEND FFTW_LIBRARIES "${FFTW_libm_LIBRARY}")
endif()


target_link_libraries(${TARGET} ${FFTW_LIBRARIES} ) 
target_link_libraries(${TARGET} m )


#
# (Optional)
# If you wish to have a "run" target for make, specify one:
#
#ADD_CUSTOM_TARGET(  run
#                    COMMAND ${NATRIUM_TARGET} 
#                    COMMENT "Run ${NATRIUM_TARGET}  with ${CMAKE_BUILD_TYPE} configuration"
#)

  
