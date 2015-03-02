
#
# define testing files
#
set(    BOLTZMANNMODELS_TEST
        boltzmannmodels/D2Q9IncompressibleModel_test.cpp
        )
set(    COLLISION_TEST
        collision/BGKTransformed_test.cpp
        collision/Collision_test.cpp
        )            
set(    PROBLEMDESCRIPTION_TEST
        problemdescription/BoundaryTools_test.cpp
        problemdescription/PeriodicBoundary2D_test.cpp
        problemdescription/MinLeeBoundary2D_test.cpp
        problemdescription/TaylorGreenTest2D.h
        problemdescription/WallTestDomain2D.h
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
		solver/SolverStats_test.cpp
		solver/ErrorStats_test.cpp
		solver/TaylorGreenVortex2D.h
		solver/TaylorGreenVortex2D.cpp
		)        
set(    TIMEINTEGRATION_TEST
        timeintegration/RungeKutta5LowStorage_test.cpp
        timeintegration/ThetaMethod_test.cpp
        timeintegration/DealIIWrapper_test.cpp
        timeintegration/ExponentialTimeIntegrator_test.cpp
        )
set(    UTILITIES_TEST
        utilities/BasicNames_test.cpp
		utilities/Logging_test.cpp
		utilities/Math_test.cpp
		utilities/CFDSolverUtilities_test.cpp
        )      



#
# define which directories to test
#
set(    TESTS
        ${BOLTZMANNMODELS_TEST}
        ${COLLISION_TEST}
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
        
set(    INTEGRATION_TEST_TARGET 
        NATriuM_Test
        )

#
# Specify a target:
#
ADD_EXECUTABLE( ${TEST_TARGET} 
                MainTest.cpp
                ${TESTS}
                )

ADD_EXECUTABLE( ${INTEGRATION_TEST_TARGET} 
                integrationtest/start-test.cpp
                integrationtest/IntegrationTestCases.h
                integrationtest/IntegrationTestCases.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-1/TaylorGreenVortex2D.h
                ${PROJECT_SOURCE_DIR}/examples/step-1/TaylorGreenVortex2D.cpp
 				${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.h
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.cpp             
                )

#
# make executable
#
deal_ii_setup_target( ${TEST_TARGET} )     
deal_ii_setup_target( ${INTEGRATION_TEST_TARGET} )
                
#
# link natrium library to executable
#           
include_directories( ${PROJECT_SOURCE_DIR}/natrium
				     ${PROJECT_SOURCE_DIR}/examples
                     ${PROJECT_BINARY_DIR}/natrium
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
                       