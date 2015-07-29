# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( intersection_2D )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if ( COMMAND cmake_policy )

cmake_policy( SET CMP0003 NEW )

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

if ( NOT CGAL_FOUND )

message(STATUS "This project requires the CGAL library, and will not be compiled.")
return()

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

message(STATUS "This project requires the Boost library, and will not be compiled.")

return()

endif()

# set directories
set ( COMMON_DIR src/Common )
set ( INTERSECTION_2D_DIR src/Intersections_2D)

# include for local directory
include_directories( ${COMMON_DIR} )

# set sources
set ( COMMON_SOURCES ${COMMON_DIR}/vector_deque.cpp)
set ( INTERSECTION_2D_SOURCES ${INTERSECTION_2D_DIR}/main.cpp ${COMMON_DIR}/intersection_functions_2.cpp ${COMMON_DIR}/triangular_mesh_2.cpp)

# Creating entries for target: intersections_2D
# ############################

add_executable( intersections_2D ${COMMON_SOURCES} ${INTERSECTION_2D_SOURCES})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_2D )

# Link the executable to CGAL and third-party libraries
target_link_libraries(intersections_2D   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

