

#
# set variable TARGET
#
set(    TARGET 
        matrix-analysis
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                analyzeMatrix.cpp
                matrixAnalysis.h
                matrixAnalysis.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-1/TaylorGreenVortex2D.h   
                ${PROJECT_SOURCE_DIR}/examples/step-1/TaylorGreenVortex2D.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.h
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.cpp
               )
               
#
# link natrium library to executable
#           
include_directories( ${PROJECT_SOURCE_DIR}/natrium 
                     ${PROJECT_BINARY_DIR}/natrium 
                     /usr/include/boost 
                   )

target_link_libraries(  ${TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
