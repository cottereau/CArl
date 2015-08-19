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
set ( SOURCE_DIR src)
set ( MESH_DIR ${SOURCE_DIR}/Triangular_mesh)
set ( INTERSECTION_DIR ${SOURCE_DIR}/Intersections)
set ( INTERSECTION_2D_EXEC_DIR ${SOURCE_DIR}/Intersections_2D_test)
set ( INTERSECTION_3D_EXEC_DIR ${SOURCE_DIR}/Intersections_3D_test)
set ( TEST_IMPORT_EXEC_DIR ${SOURCE_DIR}/Mesh_import_test)

# include for local directory
include_directories( ${SOURCE_DIR} ${MESH_DIR} ${INTERSECTION_DIR} ${INTERSECTION_2D_EXEC_DIR} ${INTERSECTION_3D_EXEC_DIR} ${MESHER_EXEC_DIR} ${TEST_IMPORT_EXEC_DIR})

# set sources
set ( INTERSECTION_2D_EXEC_SOURCES ${INTERSECTION_2D_EXEC_DIR}/main.cpp ${INTERSECTION_DIR}/intersection_functions_2.cpp ${MESH_DIR}/triangular_mesh_2.cpp)

set ( INTERSECTION_3D_EXEC_SOURCES ${INTERSECTION_3D_EXEC_DIR}/main.cpp ${INTERSECTION_DIR}/intersection_functions_3.cpp ${MESH_DIR}/triangular_mesh_3.cpp)

set ( TEST_IMPORT_EXEC_SOURCES ${TEST_IMPORT_EXEC_DIR}/main.cpp ${INTERSECTION_DIR}/intersection_functions_2.cpp ${MESH_DIR}/triangular_mesh_2.cpp)

# Creating entries for target: intersections_2D
# ############################

add_executable( intersections_2D ${INTERSECTION_2D_EXEC_SOURCES} )
add_executable( intersections_3D ${INTERSECTION_3D_EXEC_SOURCES} )
add_executable( import_mesh_2D ${TEST_IMPORT_EXEC_SOURCES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_2D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_3D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS import_mesh_2D )

set ( LIB_DIR /opt/local/lib)
set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem-mt.dylib) 

# Link the executable to CGAL and third-party libraries
target_link_libraries(intersections_2D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES} )

target_link_libraries(intersections_3D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

target_link_libraries(import_mesh_2D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
