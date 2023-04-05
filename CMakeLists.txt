cmake_minimum_required(VERSION 3.21)
project(CPPtest)

FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

set(CMAKE_CXX_STANDARD 14)

add_executable(CPPtest main.cpp)

TARGET_LINK_LIBRARIES(CPPtest LINK_PUBLIC ${Boost_LIBRARIES} )
