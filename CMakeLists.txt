cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(kth-aa)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") # -Wall -Wextra # -Werror

include_directories(include)
include_directories(src)

add_subdirectory(src)
