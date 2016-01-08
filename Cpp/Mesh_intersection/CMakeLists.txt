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

# libmesh commands
execute_process(COMMAND libmesh-config --cppflags OUTPUT_VARIABLE LIBMESH_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND libmesh-config --cxxflags OUTPUT_VARIABLE LIBMESH_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND libmesh-config --include OUTPUT_VARIABLE LIBMESH_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND libmesh-config --libs OUTPUT_VARIABLE LIBMESH_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_CPPFLAGS} ${LIBMESH_CXXFLAGS} ${LIBMESH_INCLUDE}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBMESH_INCLUDE} ${LIBMESH_LIBS}" ) 

# set directories
set ( SOURCE_DIR src )

set ( EXEC_SRC_DIR ${SOURCE_DIR}/execs )
set ( INTERSECTION_2D_EXEC_DIR ${EXEC_SRC_DIR}/intersections_2D )
set ( INTERSECTION_3D_EXEC_DIR ${EXEC_SRC_DIR}/intersections_3D )
set ( RESTRICTION_3D_EXEC_DIR  ${EXEC_SRC_DIR}/mesh_restriction_3D )
set ( TEST_IMPORT_EXEC_DIR     ${EXEC_SRC_DIR}/mesh_import_test )
set ( CARL_LIBMESH_EXEC_DIR    ${EXEC_SRC_DIR}/carl_libmesh )

set ( COMMON_SRC_DIR ${SOURCE_DIR}/common )
set ( MESH_DIR                 ${COMMON_SRC_DIR}/triangular_mesh )
set ( INTERSECTION_DIR         ${COMMON_SRC_DIR}/intersections )
set ( LIBMESH_CODE_DIR         ${COMMON_SRC_DIR}/libmesh_code )

# include for local directory
set ( LIB_DIRS ${SOURCE_DIR} ${MESH_DIR} ${INTERSECTION_DIR} ${LIBMESH_CODE_DIR} )
set ( EXEC_DIRS ${INTERSECTION_2D_EXEC_DIR} ${INTERSECTION_3D_EXEC_DIR} ${MESHER_EXEC_DIR} ${TEST_IMPORT_EXEC_DIR} ${RESTRICTION_3D_EXEC_DIR} ${CARL_LIBMESH_EXEC_DIR} )
include_directories( ${LIB_DIRS} ${EXEC_DIR} )

# set source sets
file ( GLOB MESH_2D_SOURCES "${INTERSECTION_DIR}/*2.h" "${INTERSECTION_DIR}/*2.cpp" "${MESH_DIR}/*2.cpp" "${MESH_DIR}/*2.h")
file ( GLOB MESH_3D_SOURCES "${INTERSECTION_DIR}/*3.h" "${INTERSECTION_DIR}/*3.cpp" "${MESH_DIR}/*3.cpp" "${MESH_DIR}/*3.h" "${RESTRICTION_3D_EXEC_DIR}/*3.h")
file ( GLOB LIBMESH_CODE_SOURCES "${LIBMESH_CODE_DIR}/*.h" "${LIBMESH_CODE_DIR}/*.cpp" )

# set execs sources
set ( INTERSECTION_2D_EXEC_SOURCES ${INTERSECTION_2D_EXEC_DIR}/main.cpp)
set ( INTERSECTION_3D_EXEC_SOURCES ${INTERSECTION_3D_EXEC_DIR}/main.cpp)
set ( TEST_IMPORT_EXEC_SOURCES ${TEST_IMPORT_EXEC_DIR}/main.cpp)
set ( RESTRICTION_3D_EXEC_SOURCES ${RESTRICTION_3D_EXEC_DIR}/main.cpp )
set ( CARL_LIBMESH_EXEC_SOURCES ${CARL_LIBMESH_EXEC_DIR}/main.cpp )

# Creating entries for target: intersections_2D
# ############################
add_library( MESH_3D_OBJS OBJECT ${MESH_3D_SOURCES})
add_library( MESH_2D_OBJS OBJECT ${MESH_2D_SOURCES})
add_library( LIBMESH_CODE_OBJS OBJECT ${LIBMESH_CODE_SOURCES})

add_executable( intersections_2D ${INTERSECTION_2D_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_2D_OBJS>)
add_executable( intersections_3D ${INTERSECTION_3D_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_3D_OBJS>)
add_executable( import_mesh_test ${TEST_IMPORT_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_2D_OBJS> $<TARGET_OBJECTS:MESH_3D_OBJS>)
add_executable( restrict_mesh_3D ${RESTRICTION_3D_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_3D_OBJS>)
add_executable( carl_libmesh     ${CARL_LIBMESH_EXEC_SOURCES} $<TARGET_OBJECTS:MESH_3D_OBJS> $<TARGET_OBJECTS:LIBMESH_CODE_OBJS>)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_2D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_3D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS import_mesh_test )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS restrict_mesh_3D )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS carl_libmesh )

set ( LIB_DIR /opt/local/lib)
set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem-mt.dylib) 

# Link the executable to CGAL and third-party libraries
target_link_libraries(intersections_2D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(intersections_3D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(import_mesh_test ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(restrict_mesh_3D ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(carl_libmesh ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )