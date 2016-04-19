# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Parallel_CArl )

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
message(STATUS "${CGAL_LIBRARIES}")
message(STATUS "${CGAL_VERSION}")

# libmesh commands
execute_process(COMMAND libmesh-config --cppflags OUTPUT_VARIABLE LIBMESH_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND libmesh-config --cxxflags OUTPUT_VARIABLE LIBMESH_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND libmesh-config --include OUTPUT_VARIABLE LIBMESH_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND libmesh-config --libs OUTPUT_VARIABLE LIBMESH_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_CPPFLAGS} ${LIBMESH_CXXFLAGS} ${LIBMESH_INCLUDE} -Wall" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBMESH_INCLUDE} ${LIBMESH_LIBS}" ) 

# set directories
set ( SOURCE_DIR src )

set ( EXEC_SRC_DIR ${SOURCE_DIR}/execs )
# set ( INTERSECTION_RESTRICT_3D_EXEC_DIR ${EXEC_SRC_DIR}/intersections_restrict_3D_parallel )
set ( CARL_LIBMESH_PARALLEL_EXEC_DIR    ${EXEC_SRC_DIR}/carl_libmesh_parallel )
set ( DEVEL_BACKYARD_EXEC_DIR           ${EXEC_SRC_DIR}/devel_backyard )

set ( COMMON_SRC_DIR             ${SOURCE_DIR}/common )
set ( MESH_DIR                   ${COMMON_SRC_DIR}/triangular_mesh )
set ( INTERSECTION_DIR           ${COMMON_SRC_DIR}/intersections )
set ( INTERSECTION_PARALLEL_DIR  ${COMMON_SRC_DIR}/intersections_parallel )
set ( LIBMESH_CODE_DIR           ${COMMON_SRC_DIR}/libmesh_code )

set ( ASSEMBLE_FUNCTIONS_DIR   ${EXEC_SRC_DIR}/assemble_functions )

# include for local directory
set ( LIB_DIRS ${SOURCE_DIR} ${MESH_DIR} ${INTERSECTION_DIR} ${LIBMESH_CODE_DIR} ${ASSEMBLE_FUNCTIONS_DIR} ${INTERSECTION_DIR} ${INTERSECTION_PARALLEL_DIR})
# set ( EXEC_DIRS ${INTERSECTION_RESTRICT_3D_EXEC_DIR} ${CARL_LIBMESH_PARALLEL_EXEC_DIR} ${DEVEL_BACKYARD_EXEC_DIR})
set ( EXEC_DIRS ${CARL_LIBMESH_PARALLEL_EXEC_DIR} ${DEVEL_BACKYARD_EXEC_DIR})

include_directories( ${LIB_DIRS} ${EXEC_DIR} )

# set source sets
file ( GLOB COMMON_SOURCES                "${SOURCE_DIR}/common*.h"           "${SOURCE_DIR}/common*.cpp")
file ( GLOB MESH_3D_SOURCES               "${MESH_DIR}/*3.h"                  "${MESH_DIR}/*3.cpp" "${RESTRICTION_3D_EXEC_DIR}/*3.h")
file ( GLOB INTERSECTION_SOURCES          "${INTERSECTION_DIR}/*3.h"          "${INTERSECTION_DIR}/*3.cpp")
file ( GLOB INTERSECTION_PARALLEL_SOURCES "${INTERSECTION_PARALLEL_DIR}/*.h" "${INTERSECTION_PARALLEL_DIR}/*.cpp")
file ( GLOB LIBMESH_CODE_SOURCES "${LIBMESH_CODE_DIR}/*.h" "${LIBMESH_CODE_DIR}/*.cpp" )
file ( GLOB ASSEMBLE_FUNCTIONS_SOURCES "${ASSEMBLE_FUNCTIONS_DIR}/*.h" "${ASSEMBLE_FUNCTIONS_DIR}/*.cpp" )

# set execs sources
# set ( INTERSECTION_RESTRICT_3D_EXEC_SOURCES    ${INTERSECTION_RESTRICT_3D_EXEC_DIR}/main.cpp)
set ( CARL_LIBMESH_PARALLEL_EXEC_SOURCES       ${CARL_LIBMESH_PARALLEL_EXEC_DIR}/main.cpp )
set ( DEVEL_BACKYARD_EXEC_SOURCES              ${DEVEL_BACKYARD_EXEC_DIR}/main.cpp )

# Creating entries for target: intersections_2D
# ############################
add_library( COMMON_OBJS OBJECT                         ${COMMON_SOURCES})
add_library( MESH_3D_OBJS OBJECT                        ${MESH_3D_SOURCES})
add_library( INTERSECTION_OBJS OBJECT                   ${INTERSECTION_SOURCES})
add_library( INTERSECTION_PARALLEL_OBJS OBJECT          ${INTERSECTION_PARALLEL_SOURCES})
add_library( LIBMESH_CODE_OBJS OBJECT                   ${LIBMESH_CODE_SOURCES})
add_library( ASSEMBLE_FUNCTIONS_OBJS OBJECT             ${ASSEMBLE_FUNCTIONS_SOURCES})

# add_executable( intersections_restrict_3D_parallel ${INTERSECTION_RESTRICT_3D_EXEC_SOURCES} 
#                                                    $<TARGET_OBJECTS:MESH_3D_OBJS> 
#                                                    $<TARGET_OBJECTS:COMMON_OBJS>)

add_executable( carl_libmesh_parallel     ${CARL_LIBMESH_PARALLEL_EXEC_SOURCES} 
                                                  $<TARGET_OBJECTS:MESH_3D_OBJS> 
                                                  $<TARGET_OBJECTS:INTERSECTION_OBJS> 
                                                  $<TARGET_OBJECTS:LIBMESH_CODE_OBJS> 
                                                  $<TARGET_OBJECTS:COMMON_OBJS>
                                                  $<TARGET_OBJECTS:ASSEMBLE_FUNCTIONS_OBJS>)
                                                  
add_executable( test_devel_parallel     ${DEVEL_BACKYARD_EXEC_SOURCES} 
                                                  $<TARGET_OBJECTS:MESH_3D_OBJS> 
                                                  $<TARGET_OBJECTS:INTERSECTION_PARALLEL_OBJS> 
                                                  $<TARGET_OBJECTS:LIBMESH_CODE_OBJS> 
                                                  $<TARGET_OBJECTS:COMMON_OBJS>
                                                  $<TARGET_OBJECTS:ASSEMBLE_FUNCTIONS_OBJS>)
                                                  
# add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_restrict_3D_parallel )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS carl_libmesh_parallel )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS test_libmesh_parallel )

set ( LIB_DIR /opt/local/lib)
set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem-mt.dylib) 

# Link the executable to CGAL and third-party libraries
# target_link_libraries(intersections_restrict_3D_parallel ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(carl_libmesh_parallel ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(test_devel_parallel  ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )