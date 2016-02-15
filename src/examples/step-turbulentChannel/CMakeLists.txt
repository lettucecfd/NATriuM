

#
# set variable TARGET
#
set(    TARGET 
        step-turbulent-channel
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                step-turbulent-channel.cpp
                TurbulentChannelFlow3D.h   
                TurbulentChannelFlow3D.cpp
               )
#
# link natrium library to executable
#           
include_directories( /usr/include/boost ${PROJECT_SOURCE_DIR}/library) 

target_link_libraries(  ${TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
