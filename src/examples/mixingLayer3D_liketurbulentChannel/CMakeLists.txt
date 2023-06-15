# set variable TARGET
#
set(TARGET mixingLayer_likeChannel)
# Specify a target:
ADD_EXECUTABLE( ${TARGET}
        step-mixingLayer.cpp
        MixingLayer3D.h
        MixingLayer3D.cpp
        FinalChannelStatistics.h
        FinalChannelStatistics.cpp
        )
# link natrium library to executable
include_directories( /usr/include/boost ${PROJECT_SOURCE_DIR}/library)
target_link_libraries(  ${TARGET}
        natrium
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        )
# make executable
deal_ii_setup_target( ${TARGET} )