

#
# set variable TARGET
#
set(    TARGET 
	shear-layer    
	#shear-layer-MinionBrown
        )
 

#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                shear_layer.cpp
		#                shear_layer_MinionBrown.cpp
                EnstrophySubdomain.h
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
        

                
