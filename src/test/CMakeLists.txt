#### include(../CMake_UnityBuild.txt)


file(   GLOB_RECURSE
        TEST_FILES
        advection/*.cpp
        boundaries/*.cpp
        collision/*.cpp
        problemdescription/*.cpp
        smoothing/*.cpp
        solver/*.cpp
        stencils/*.cpp
        timeintegration/*.cpp
        utilities/*.cpp
        )
      
      
#
# set variable TEST_TARGET (only used in this file)
#
set(    TEST_TARGET 
        NATriuM_UnitTest_exe
        )
        
set(    INTEGRATION_TEST_TARGET 
        NATriuM_Test
        ) 
 


#enable_unity_build(${TEST_TARGET} TEST_FILES)

#
# Specify a target:
#
ADD_EXECUTABLE( ${TEST_TARGET} 
                MainTest.cpp
                ${TEST_FILES}
                )

ADD_EXECUTABLE( ${INTEGRATION_TEST_TARGET} 
                integrationtest/start-test.cpp
                integrationtest/IntegrationTestCases.h
                integrationtest/IntegrationTestCases.cpp       
                )

#
# make executable
#
deal_ii_setup_target( ${TEST_TARGET} )     
deal_ii_setup_target( ${INTEGRATION_TEST_TARGET} )
                
#
# link natrium library to executable
#           
include_directories( ${PROJECT_SOURCE_DIR}/library
                     /usr/include/boost
                     )

target_link_libraries(  ${TEST_TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                        ) 
                        
target_link_libraries(  ${INTEGRATION_TEST_TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                        )                        
                       
