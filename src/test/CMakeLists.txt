
#
# define testing files
#
set(    BOLTZMANNMODELS_TEST
        boltzmannmodels/D2Q9IncompressibleModel_test.cpp
        )
set(    COLLISIONMODELS_TEST
        collisionmodels/BGKTransformed_test.cpp
        )            
set(    PROBLEMDESCRIPTION_TEST
        problemdescription/PeriodicBoundary1D_test.cpp
        problemdescription/BoundaryTools_test.cpp
        )   
set(    ADVECTION_TEST
        advection/SEDGMinLee_test.cpp
        advection/PeriodicTestDomain2D.h
        advection/PeriodicTestDomain2D.cpp
        )     
set(	SOLVER_TEST
		solver/DistributionFunctions_test.cpp
		solver/CFDSolver_test.cpp
		solver/PeriodicTestFlow2D.h
		solver/SolverConfiguration_test.cpp	
		)        
set(    TIMEINTEGRATION_TEST
        timeintegration/RungeKutta5LowStorage_test.cpp
        )
set(    UTILITIES_TEST
        utilities/BasicNames_test.cpp
        )      



#
# define which directories to test
#
set(    TESTS
        ${BOLTZMANNMODELS_TEST}
        ${COLLISIONMODELS_TEST}
        ${PROBLEMDESCRIPTION_TEST}
        ${ADVECTION_TEST}
        ${SOLVER_TEST}
        ${TIMEINTEGRATION_TEST}
        ${UTILITIES_TEST}
                
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
include_directories( ${PROJECT_SOURCE_DIR}/natrium
                     ${PROJECT_BINARY_DIR}/natrium
                     /usr/include/boost
                     )

target_link_libraries(  ${TEST_TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                        )
                       