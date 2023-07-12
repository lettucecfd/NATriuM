# set variable TARGET
set(TARGET step-mixingLayer_solenoidal)
# Specify a target:
ADD_EXECUTABLE( ${TARGET}
        step-mixingLayer_solenoidal.cpp
        MixingLayer3D_solenoidal.h
        MixingLayer3D_solenoidal.cpp
        )
# link natrium library to executable
include_directories( /usr/include/boost ../../../library)
target_link_libraries(  ${TARGET}
        natrium
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        )
# make executable
deal_ii_setup_target( ${TARGET} )