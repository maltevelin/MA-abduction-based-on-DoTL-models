cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 20)

project(abduction)

include_directories( ${PROJECT_SOURCE_DIR}/src )

file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*")

set(CMAKE_BUILD_TYPE Debug)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

add_executable( ${PROJECT_NAME} ${SOURCES} )
