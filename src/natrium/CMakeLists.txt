
 
#
# Store all library files in environment variable SRC_FILES
#

file(   GLOB_RECURSE
        SRC_FILES
        benchmarks/*
        boltzmannmodels/*
        collisionmodels/*
        problemdescription/*
        solver/*
        streamingdata/*
        timeintegration/*
        utilites/*
        )
 
 
 
#
# Create natrium library
#
add_library( natrium ${SRC_FILES} )



#
# Make library
#
DEAL_II_SETUP_TARGET ( natrium )


#
# Include src directory (to avoid ../../)
#
include_directories( ${PROJECT_SOURCE_DIR}/natrium )

#
# copy natrium library to lib/ folder when typing 'make install'
#
install (TARGETS natrium DESTINATION ${PROJECT_BINARY_DIR}/lib)
                        
               