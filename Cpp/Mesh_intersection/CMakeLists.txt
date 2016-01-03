# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( intersection_2D )

cmake_minimum_required(VERSION 2.8.8)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.8)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.8)
    cmake_policy(VERSION 2.8.8)
  endif()
endif()


set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if ( COMMAND cmake_policy )

cmake_policy( SET CMP0003 NEW )

endif()

# CGAL and its components
find_package( CGAL COMPONENTS Core )

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

message(STATUS "${CMAKE_MODULE_PATH}")

# set directories
set ( SOURCE_DIR src )

set ( EXEC_SRC_DIR ${SOURCE_DIR}/execs )
set ( INTERSECTION_2D_EXEC_DIR ${EXEC_SRC_DIR}/intersections_2D )
set ( INTERSECTION_3D_EXEC_DIR ${EXEC_SRC_DIR}/intersections_3D )
set ( RESTRICTION_3D_EXEC_DIR  ${EXEC_SRC_DIR}/mesh_restriction_3D )
set ( TEST_IMPORT_EXEC_DIR     ${EXEC_SRC_DIR}/mesh_import_test )

set ( COMMON_SRC_DIR ${SOURCE_DIR}/common )
set ( MESH_DIR                 ${COMMON_SRC_DIR}/triangular_mesh )
set ( INTERSECTION_DIR         ${COMMON_SRC_DIR}/intersections )

# include for local directory
include_directories( ${SOURCE_DIR} ${MESH_DIR} ${INTERSECTION_DIR} ${INTERSECTION_2D_EXEC_DIR} ${INTERSECTION_3D_EXEC_DIR} ${MESHER_EXEC_DIR} ${TEST_IMPORT_EXEC_DIR} ${RESTRICTION_3D_EXEC_DIR} )

# set source sets
set ( MESH_2D_SOURCES ${INTERSECTION_DIR}/intersection_functions_2.cpp ${INTERSECTION_DIR}/triangular_mesh_intersection_2.cpp ${MESH_DIR}/triangular_mesh_2.cpp )
set ( MESH_3D_SOURCES ${INTERSECTION_DIR}/intersection_functions_3.cpp ${INTERSECTION_DIR}/triangular_mesh_intersection_3.cpp ${MESH_DIR}/triangular_mesh_3.cpp )

# set execs sources
set ( INTERSECTION_2D_EXEC_SOURCES ${INTERSECTION_2D_EXEC_DIR}/main.cpp)
set ( INTERSECTION_3D_EXEC_SOURCES ${INTERSECTION_3D_EXEC_DIR}/main.cpp)
set ( TEST_IMPORT_EXEC_SOURCES ${TEST_IMPORT_EXEC_DIR}/main.cpp)
set ( RESTRICTION_3D_EXEC_SOURCES ${RESTRICTION_3D_EXEC_DIR}/main.cpp ${RESTRICTION_3D_EXEC_DIR}/mesh_restriction_3D.h )

# Creating entries for target: intersections_2D
# ############################
add_library( MESH_3D_OBJS OBJECT ${MESH_3D_SOURCES})
add_library( MESH_2D_OBJS OBJECT ${MESH_2D_SOURCES})

add_executable( intersections_2D ${INTERSECTION_2D_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_2D_OBJS>)
add_executable( intersections_3D ${INTERSECTION_3D_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_3D_OBJS>)
add_executable( import_mesh_test ${TEST_IMPORT_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_2D_OBJS> $<TARGET_OBJECTS:MESH_3D_OBJS>)
add_executable( restrict_mesh_3D ${RESTRICTION_3D_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_3D_OBJS>)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_2D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_3D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS import_mesh_test )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS restrict_mesh_3D )

set ( LIB_DIR /opt/local/lib)
set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem-mt.dylib) 

# Link the executable to CGAL and third-party libraries
target_link_libraries(intersections_2D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(intersections_3D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(import_mesh_test ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(restrict_mesh_3D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
