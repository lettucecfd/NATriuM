

#
# set variable TARGET
#
set(    TARGET 
        convergence-analysis-wall-p
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                convergence-analysis-wall-p.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.h   
                ${PROJECT_SOURCE_DIR}/examples/step-2/CouetteFlow2D.cpp
               )
#
# link natrium library to executable
#           
include_directories( /usr/include/boost  ${PROJECT_SOURCE_DIR}/library)

target_link_libraries(  ${TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
