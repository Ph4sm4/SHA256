cmake_minimum_required(VERSION 3.11.0)
project(sha256 VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")

set(PROJECT_SOURCES main.cpp sha256.h)

# Define the main executable
add_executable(sha256 ${PROJECT_SOURCES})