

#
# set variable TARGET
#
set(    TARGET 
        step-3-LidDrivenCavity2D
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                step-3.cpp
                LidDrivenCavity2D.h   
                LidDrivenCavity2D.cpp
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
        

                
