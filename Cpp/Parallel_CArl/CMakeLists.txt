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

message(STATUS "CGAL version : ${CGAL_VERSION}")
message(STATUS "CGAL dir     : ${CGAL_DIR}")

# libmesh commands
if ( NOT LIBMESH_DIR )
   set (LIBMESH_DIR "/usr/local/lib")
endif()

if ( NOT LIBMESH_CONFIG_DIR )
   set (LIBMESH_CONFIG_DIR "/usr/local/bin")
endif()

message (STATUS "Using ${LIBMESH_DIR} as the libMesh directory")

set ( LIBMESH_CONFIG ${LIBMESH_CONFIG_DIR}/libmesh-config )
execute_process(COMMAND ${LIBMESH_CONFIG} --cppflags OUTPUT_VARIABLE LIBMESH_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --cxxflags OUTPUT_VARIABLE LIBMESH_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --include OUTPUT_VARIABLE LIBMESH_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --libs OUTPUT_VARIABLE LIBMESH_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_CPPFLAGS} ${LIBMESH_CXXFLAGS} ${LIBMESH_INCLUDE} -Wall" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBMESH_INCLUDE} ${LIBMESH_LIBS}" ) 

# set directories
set ( SOURCE_DIR src )

set ( EXEC_SRC_DIR ${SOURCE_DIR}/execs )
set ( CARL_LIBMESH_PARALLEL_EXEC_DIR    ${EXEC_SRC_DIR}/carl_libmesh_parallel )
set ( DEVEL_BACKYARD_EXEC_DIR           ${EXEC_SRC_DIR}/devel_backyard )
set ( TEST_CASES_EXEC_DIR               ${EXEC_SRC_DIR}/tests )
set ( ASSEMBLE_FUNCTIONS_DIR            ${EXEC_SRC_DIR}/assemble_functions )

set ( COMMON_SRC_DIR             ${SOURCE_DIR}/common )

set ( MESH_DIR                   ${COMMON_SRC_DIR}/triangular_mesh )
set ( INTERSECTION_PARALLEL_DIR  ${COMMON_SRC_DIR}/intersections_parallel )
set ( LIBMESH_CODE_DIR           ${COMMON_SRC_DIR}/libmesh_code )

set ( LM_BASE_CG_SOLVER_DIR     ${LIBMESH_CODE_DIR}/base_CG_solver )
set ( LM_COUPLED_SOLVER_DIR     ${LIBMESH_CODE_DIR}/coupled_solver )
set ( LM_COUPLED_SYSTEM_DIR     ${LIBMESH_CODE_DIR}/coupled_system )
set ( LM_PETSC_MATRIX_OPS_DIR   ${LIBMESH_CODE_DIR}/PETSC_matrix_operations )
set ( LM_SYSTEM_SOLVER_DIR      ${LIBMESH_CODE_DIR}/system_solver )
set ( LM_MISC_DIR               ${LIBMESH_CODE_DIR}/misc )

# include for local directory
set ( LIBMESH_LIB_DIRS ${LM_BASE_CG_SOLVER_DIR} ${LM_COUPLED_SOLVER_DIR} ${LM_COUPLED_SYSTEM_DIR} ${LM_PETSC_MATRIX_OPS_DIR} ${LM_SYSTEM_SOLVER_DIR} ${LM_MISC_DIR})
set ( LIB_DIRS ${SOURCE_DIR} ${MESH_DIR} ${INTERSECTION_DIR} ${ASSEMBLE_FUNCTIONS_DIR} ${INTERSECTION_PARALLEL_DIR})
set ( EXEC_DIRS ${CARL_LIBMESH_PARALLEL_EXEC_DIR} ${DEVEL_BACKYARD_EXEC_DIR} ${TEST_CASES_EXEC_DIR})

include_directories( ${LIBMESH_LIB_DIRS} ${LIB_DIRS} ${EXEC_DIR} )

# set source sets
file ( GLOB COMMON_SOURCES                "${COMMON_SRC_DIR}/common*.h"           
                                          "${COMMON_SRC_DIR}/common*.cpp")
file ( GLOB MESH_3D_SOURCES               "${MESH_DIR}/*3.h"                  
                                          "${MESH_DIR}/*3.cpp"
                                          "${RESTRICTION_3D_EXEC_DIR}/*3.h")
file ( GLOB INTERSECTION_PARALLEL_SOURCES "${INTERSECTION_PARALLEL_DIR}/*.h"
                                          "${INTERSECTION_PARALLEL_DIR}/*.cpp")
file ( GLOB ASSEMBLE_FUNCTIONS_SOURCES    "${ASSEMBLE_FUNCTIONS_DIR}/*.h"
                                          "${ASSEMBLE_FUNCTIONS_DIR}/*.cpp" )

file ( GLOB LM_BASE_CG_SOLVER_SOURCES     "${LM_BASE_CG_SOLVER_DIR}/*.h"
                                          "${LM_BASE_CG_SOLVER_DIR}/*.cpp" )
file ( GLOB LM_COUPLED_SOLVER_SOURCES     "${LM_COUPLED_SOLVER_DIR}/*.h"
                                          "${LM_COUPLED_SOLVER_DIR}/*.cpp" )
file ( GLOB LM_COUPLED_SYSTEM_SOURCES     "${LM_COUPLED_SYSTEM_DIR}/*.h"
                                          "${LM_COUPLED_SYSTEM_DIR}/*.cpp" )
file ( GLOB LM_PETSC_MATRIX_OPS_SOURCES   "${LM_PETSC_MATRIX_OPS_DIR}/*.h"
                                          "${LM_PETSC_MATRIX_OPS_DIR}/*.cpp" )
file ( GLOB LM_SYSTEM_SOLVER_SOURCES      "${LM_SYSTEM_SOLVER_DIR}/*.h"
                                          "${LM_SYSTEM_SOLVER_DIR}/*.cpp" )
file ( GLOB LM_MISC_SOURCES               "${LM_MISC_DIR}/*.h"
                                          "${LM_MISC_DIR}/*.cpp" )

# set execs sources
set ( CARL_LIBMESH_PARALLEL_EXEC_SOURCES       ${CARL_LIBMESH_PARALLEL_EXEC_DIR}/main.cpp )
set ( DEVEL_BACKYARD_EXEC_SOURCES              ${DEVEL_BACKYARD_EXEC_DIR}/main.cpp )
file ( GLOB TEST_EXECS ${TEST_CASES_EXEC_DIR}/*.cpp )

# Creating entries for target: intersections_2D
# ############################
add_library( COMMON_OBJS OBJECT                         "${SOURCE_DIR}/carl_headers.h" "${COMMON_SRC_DIR}/CGAL_typedefs.h" ${COMMON_SOURCES})
add_library( MESH_3D_OBJS OBJECT                        ${MESH_3D_SOURCES})
add_library( INTERSECTION_PARALLEL_OBJS OBJECT          ${INTERSECTION_PARALLEL_SOURCES})
add_library( LIBMESH_CODE_OBJS OBJECT                   ${LM_BASE_CG_SOLVER_SOURCES} ${LM_COUPLED_SOLVER_SOURCES} ${LM_COUPLED_SYSTEM_SOURCES} ${LM_PETSC_MATRIX_OPS_SOURCES} ${LM_SYSTEM_SOLVER_SOURCES} ${LM_MISC_SOURCES})
add_library( ASSEMBLE_FUNCTIONS_OBJS OBJECT             ${ASSEMBLE_FUNCTIONS_SOURCES})

add_executable( carl_libmesh_parallel     ${CARL_LIBMESH_PARALLEL_EXEC_SOURCES} 
                                                  $<TARGET_OBJECTS:MESH_3D_OBJS> 
                                                  $<TARGET_OBJECTS:LIBMESH_CODE_OBJS> 
                                                  $<TARGET_OBJECTS:COMMON_OBJS>
                                                  $<TARGET_OBJECTS:ASSEMBLE_FUNCTIONS_OBJS>)
                                                  
add_executable( test_devel_parallel     ${DEVEL_BACKYARD_EXEC_SOURCES} 
                                                  $<TARGET_OBJECTS:MESH_3D_OBJS> 
                                                  $<TARGET_OBJECTS:INTERSECTION_PARALLEL_OBJS> 
                                                  $<TARGET_OBJECTS:LIBMESH_CODE_OBJS> 
                                                  $<TARGET_OBJECTS:COMMON_OBJS>
                                                  $<TARGET_OBJECTS:ASSEMBLE_FUNCTIONS_OBJS>)

# Set up the test files                                           
foreach( test_source_file ${TEST_EXECS} )
    string ( REPLACE ".cpp" "" full_test_exec_name ${test_source_file} )
    string ( REPLACE "${CMAKE_CURRENT_LIST_DIR}/${TEST_CASES_EXEC_DIR}/" "" test_exec_name ${full_test_exec_name} )
    add_executable ( ${test_exec_name} ${test_source_file} $<TARGET_OBJECTS:MESH_3D_OBJS> 
                                                           $<TARGET_OBJECTS:INTERSECTION_PARALLEL_OBJS> 
                                                           $<TARGET_OBJECTS:LIBMESH_CODE_OBJS> 
                                                           $<TARGET_OBJECTS:COMMON_OBJS>
                                                           $<TARGET_OBJECTS:ASSEMBLE_FUNCTIONS_OBJS>)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${test_exec_name} )
    target_link_libraries( ${test_exec_name} ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
endforeach( test_source_file ${TEST_EXECS} )

                                                  
# add_to_cached_list( CGAL_EXECUTABLE_TARGETS intersections_restrict_3D_parallel )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS carl_libmesh_parallel )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS test_libmesh_parallel )

# Link the executable to CGAL and third-party libraries
target_link_libraries(carl_libmesh_parallel ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(test_devel_parallel  ${CGAL_LIBRARIES} ${BOOST_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )