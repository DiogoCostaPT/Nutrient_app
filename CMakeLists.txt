cmake_minimum_required(VERSION 3.10)

IF (UNIX AND NOT APPLE)
	MESSAGE(STATUS "current platform : Linux")
ELSEIF(APPLE)	
	MESSAGE(STATUS "current platform : MacOS")
ELSE()
	MESSAGE(STATUS "current platform : Windows")
ENDIF()
	
file(MAKE_DIRECTORY dist/Debug/GNU-Linux)
	
project(nutrient_app)
file(GLOB SOURCES "main.c")
	
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-g -Wall -lm")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/dist/Debug/GNU-Linux)
	
LINK_LIBRARIES(m)
add_executable(nutrient_app ${SOURCES})

	
