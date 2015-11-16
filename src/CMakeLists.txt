#
# CMake File for NATriuM test
#    

    
#
# It is good practise to specify a version requirement:
#
cmake_minimum_required( VERSION 2.8.8 )

#
# Find deal library
#
find_package(   deal.II 8.0 REQUIRED
                HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)

#
# Initialize Deal variables (has to be called before "project")
#
DEAL_II_INITIALIZE_CACHED_VARIABLES() 
 
#
# PROJECT has to be called after DEAL_II_INITIALIZE_CACHED_VARIABLES
# otherwise CMake will have already set the compiler and compiler
# flags.
#
project( NATriuM )
 
#
# include boost
#
SET(BOOST_MIN_VERSION "1.54.0")
find_package( Boost ${BOOST_MIN_VERSION}
              COMPONENTS
              system
              filesystem
              unit_test_framework
              REQUIRED
              )
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
              
#
# descend into subdirectories (build natrium library in src and unit tests in test as well as examples)
#
add_subdirectory( library )
add_subdirectory( test )
#add_subdirectory( examples )
add_subdirectory( analysis )
 
#
# set variable NATRIUM_TARGET (only used in this file)
#
set(    NATRIUM_TARGET 
        #${PROJECT_BINARY_DIR}/NATriuM_exe
        NATriuM_exe
        )


  
#
# Specify a target:
#
add_executable( ${NATRIUM_TARGET} 
                NATriuM.cpp
                )

#
# make executable
#
deal_ii_setup_target( ${NATRIUM_TARGET} )     
         
                
#
# link natrium library to executable
#           
include_directories( "${PROJECT_SOURCE_DIR}/library" )
target_link_libraries(  ${NATRIUM_TARGET}
                        natrium
                        )
                       
 

#
# (Optional)
# If you wish to have a "run" target for make, specify one:
#
ADD_CUSTOM_TARGET(  run
                    COMMAND ${NATRIUM_TARGET} 
                    COMMENT "Run ${NATRIUM_TARGET}  with ${CMAKE_BUILD_TYPE} configuration"
)
ADD_CUSTOM_TARGET(  test
                    COMMAND ./bin/test/NATriuM_UnitTest_exe 
                    COMMENT "Run UnitTests  with ${CMAKE_BUILD_TYPE} configuration"
)


  
