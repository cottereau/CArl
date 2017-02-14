project( Libmesh_test )

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

# Boost and its components
find_package( Boost COMPONENTS serialization REQUIRED )

if ( NOT Boost_FOUND )
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()

SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ansi -pedantic -O3 -std=gnu++11" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_compact_unwind" )

message(STATUS "${CMAKE_CXX_FLAGS}")

# include for Boost directory
include_directories(${Boost_INCLUDE_DIRS})

# includes for Voro++
set ( VORO_PP_HEADERS "/Users/breubreubreu/Execs/voro++-0.4.6/voro/include/voro++" )
set ( VORO_PP_LIBS "/Users/breubreubreu/Execs/voro++-0.4.6/voro/lib" )

include_directories ( ${VORO_PP_HEADERS} )
link_directories ( ${VORO_PP_LIBS} )

# set directories
set ( SOURCE_DIR src)

# include for local directory
include_directories( ${SOURCE_DIR} )

# set source sets
set ( VORO_LIBS "${SOURCE_DIR}/voro_functions.h" "${SOURCE_DIR}/voro_functions.cpp" )
set ( VORO_PARSER_LIBS "${SOURCE_DIR}/raw_data_parser.h" "${SOURCE_DIR}/raw_data_parser.cpp" )
set ( VORO_TEST_SRC "${SOURCE_DIR}/getpot.h" "${SOURCE_DIR}/main.cpp" ) 
set ( VORO_TEST_CUT_SRC "${SOURCE_DIR}/getpot.h" "${SOURCE_DIR}/test_semi_plan_wall.cpp" ) 
set ( VORO_TEST_SEMI_SRC "${SOURCE_DIR}/getpot.h" "${SOURCE_DIR}/main_walls.cpp" ) 
set ( VORO_PARSER "${SOURCE_DIR}/getpot.h" "${SOURCE_DIR}/parse_test.cpp" )
# Creating entries
# ############################

add_executable( voro3D ${VORO_TEST_SRC} ${VORO_LIBS} )
add_executable( voro3D_half ${VORO_TEST_SEMI_SRC} ${VORO_LIBS} )
add_executable( test_half ${VORO_TEST_CUT_SRC} )
add_executable (data_parser ${VORO_PARSER} ${VORO_LIBS} ${VORO_PARSER_LIBS})

# set ( LIB_DIR /opt/local/lib)
# set ( BOOST_LIBRARIES ${LIB_DIR}/libboost_filesystem-mt.dylib ${LIB_DIR}/libboost_serialization-mt.dylib)

# Link the executable
target_link_libraries( voro3D ${BOOST_LIBRARIES} ${Boost_LIBRARIES} voro++ )
target_link_libraries( voro3D_half ${BOOST_LIBRARIES} ${Boost_LIBRARIES} voro++ )
target_link_libraries( test_half ${BOOST_LIBRARIES} ${Boost_LIBRARIES} voro++ )
target_link_libraries( data_parser ${BOOST_LIBRARIES} ${Boost_LIBRARIES} voro++ )
