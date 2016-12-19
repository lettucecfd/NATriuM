

#
# set variable TARGET
#
set(    TARGET 
        step-16-taylorGreenVortex3D
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                step-16.cpp
               )
#
# link natrium library to executable
#           
include_directories( /usr/include/boost ${PROJECT_SOURCE_DIR}/library )

target_link_libraries(  ${TARGET}
                        natrium
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_PROGRAM_OPTIONS_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
