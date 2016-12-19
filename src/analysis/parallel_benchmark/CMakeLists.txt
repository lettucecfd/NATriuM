

#
# set variable TARGET
#
set(    TARGET 
        benchmark    
        )
 

#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                benchmark.cpp
               )
               
#
# link natrium library to executable
#           
include_directories( /usr/include/boost  ${PROJECT_SOURCE_DIR}/library) 

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
        

                
