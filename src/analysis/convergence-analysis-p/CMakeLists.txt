

#
# set variable TARGET
#
set(    TARGET 
        convergence-analysis-p
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                convergence-analysis-p.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-1/TaylorGreenVortex2D.h   
                ${PROJECT_SOURCE_DIR}/examples/step-1/TaylorGreenVortex2D.cpp
               )
#
# link natrium library to executable
#           
include_directories( /usr/include/boost ${PROJECT_SOURCE_DIR}/natrium ${PROJECT_BINARY_DIR}/natrium)

target_link_libraries(  ${TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
