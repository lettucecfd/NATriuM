

#
# set variable TARGET
#
set(    TARGET 
        matrix-analysis
        )


#
# Specify a target:
#
ADD_EXECUTABLE( ${TARGET} 
                analyzeMatrix.cpp
                matrixAnalysis.h
                matrixAnalysis.cpp
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
                       ${Boost_PROGRAM_OPTIONS_LIBRARY}
                        )
                       
                       
#
# make executable
#
deal_ii_setup_target( ${TARGET} )     
        

                
