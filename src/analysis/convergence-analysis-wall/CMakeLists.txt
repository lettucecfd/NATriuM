

#
# set variable TARGET
#
set(    TARGET 
        convergence-analysis-wall
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                convergence-analysis-wall.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.h   
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.cpp
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
        

                
