

#
# set variable TARGET
#
set(    TARGET 
        step-1-couette2D     
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                step-1.cpp
                CouetteFlow2D.h   
                CouetteFlow2D.cpp
               )
#
# link natrium library to executable
#           
include_directories( /usr/include/boost ../../bin/src ../../src ${PROJECT_BINARY_DIR}/src)

target_link_libraries(  ${TARGET}
                        natrium
                      #${Boost_FILESYSTEM_LIBRARY}
                      # ${Boost_SYSTEM_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
