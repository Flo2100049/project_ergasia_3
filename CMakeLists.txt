# Minimum required version of CMake
cmake_minimum_required(VERSION 3.1...3.23)

# Project name
project(main)

# Find CGAL and its components
find_package(CGAL QUIET COMPONENTS Qt5)
find_package(Boost 1.75 REQUIRED COMPONENTS json)

if (NOT CGAL_FOUND)
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# Find Boost
find_package(Boost REQUIRED)

if (NOT Boost_FOUND)
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Include the directory where custom headers are located
#target_include_directories(main PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Add the include directory for custom headers
include_directories(include)

# Create the executable
add_executable(main main.cpp)

set(CMAKE_BUILD_TYPE Release)

# Link the executable to CGAL and Qt5 libraries
target_link_libraries(main PRIVATE CGAL::CGAL)  
target_link_libraries(main PRIVATE Boost::json)

if(CGAL_Qt5_FOUND)
   add_definitions(-DCGAL_USE_BASIC_VIEWER)
   target_link_libraries(main PRIVATE CGAL::CGAL CGAL::CGAL_Qt5)
endif()