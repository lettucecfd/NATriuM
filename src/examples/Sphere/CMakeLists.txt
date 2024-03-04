

#
# set variable TARGET
#
set(    TARGET 
        step-sphere
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET}
        step-sphere.cpp
        Sphere.h
        Sphere.cpp
               )
#
# link natrium library to executable
#           
include_directories( /usr/include/boost ../../library)

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
        

                
