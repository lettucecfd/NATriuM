

#
# set variable TARGET
#
set(    TARGET 
        sinus-shear-analysis
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                sinus-shear.cpp
                ${PROJECT_SOURCE_DIR}/examples/step-sinusoidal/SinusoidalShear2D.h   
                ${PROJECT_SOURCE_DIR}/examples/step-sinusoidal/SinusoidalShear2D.cpp             
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
        

                
