
#
# define testing files
#
set(    BOLTZMANNMODELS_TEST
        boltzmannmodels/D2Q9IncompressibleModel_test.cpp
        )
set(    COLLISIONMODELS_TEST
        collisionmodels/BGKTransformed_test.cpp
        )        
set(    BENCHMARKS_TEST
        benchmarks/CouetteFlow2D_test.cpp
        )        



#
# define which directories to test
#
set(    TESTS
        ${BENCHMARKS_TEST}
        ${BOLTZMANNMODELS_TEST}
        ${COLLISIONMODELS_TEST}
        
        
)


#
# set variable TEST_TARGET (only used in this file)
#
set(    TEST_TARGET 
        NATriuM_UnitTest_exe
        )

#
# Specify a target:
#
ADD_EXECUTABLE( ${TEST_TARGET} 
                MainTest.cpp
                ${TESTS}
                )

#
# make executable
#
deal_ii_setup_target( ${TEST_TARGET} )     
         
                
#
# link natrium library to executable
#           
include_directories( ${PROJECT_SOURCE_DIR}/src
                     ${PROJECT_BINARY_DIR}/src
                     /usr/include/boost
                     )

target_link_libraries(  ${TEST_TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                        )
                       