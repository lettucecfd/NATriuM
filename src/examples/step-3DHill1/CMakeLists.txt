

#
# set variable TARGET
#
set(    TARGET 
        step-hill3D
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                step-hill3D.cpp
                Hill3D.h   
                Hill3D.cpp
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
        

                
