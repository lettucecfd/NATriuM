

#
# set variable TARGET
#
set(    TARGET 
        step-turbulent-channel
        )

set (	RESTART_TARGET
		convert-old-turbulent-channel-checkpoint
		)

# 
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                step-turbulent-channel.cpp
                TurbulentChannelFlow3D.h   
                TurbulentChannelFlow3D.cpp
                FinalChannelStatistics.h
                FinalChannelStatistics.cpp
               )
               
ADD_EXECUTABLE( ${RESTART_TARGET} 
                convert-old-turbulent-channel-checkpoint.cpp
                TurbulentChannelFlow3D.h   
                TurbulentChannelFlow3D.cpp
                FinalChannelStatistics.h
                FinalChannelStatistics.cpp
               )
#

# link natrium library to executable
#           
include_directories( /usr/include/boost ${PROJECT_SOURCE_DIR}/library) 

target_link_libraries(  ${TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_PROGRAM_OPTIONS_LIBRARY}
                        )
target_link_libraries(  ${RESTART_TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_PROGRAM_OPTIONS_LIBRARY}
                        )      
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
deal_ii_setup_target( ${RESTART_TARGET} )     

                
