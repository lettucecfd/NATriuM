include(../CMake_UnityBuild.txt)

 
#
# Store all library files in environment variable SRC_FILES
#

file( GLOB_RECURSE 
      INCLUDES
      natrium/*.h
      )

file(   GLOB_RECURSE
        SRC_FILES
        natrium/advection/*.cpp
        natrium/collision/*.cpp
        natrium/problemdescription/*.cpp
        natrium/solver/*.cpp
        natrium/stencils/*.cpp
        natrium/timeintegration/*.cpp
        natrium/utilities/*.cpp
        )


        
enable_unity_build(natrium SRC_FILES)
  
 
 
#
# Create natrium library
#
add_library( natrium SHARED ${INCLUDES} ${SRC_FILES} )



#
# Make library
#
DEAL_II_SETUP_TARGET ( natrium )


#
# Include src directory (to avoid ../../)
#
# include_directories( ${PROJECT_SOURCE_DIR}/natrium )

#
# copy natrium library to lib/ folder when typing 'make install'
#
install (TARGETS natrium DESTINATION ${PROJECT_BINARY_DIR}/lib)
                        
               