

#
# set variable TARGET
#
set(    TARGET 
		sl-stability        
    )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                sl-stability.cpp
               )
               
#
# link natrium library to executable
#           
include_directories(  ${PROJECT_SOURCE_DIR}/library
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
        

                
