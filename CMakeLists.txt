cmake_minimum_required(VERSION 3.6)
project(ExtendedKF)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/FusionEKF.cpp
        src/kalman_filter.cpp
        src/main.cpp
        src/tools.cpp)

add_executable(ExtendedKF ${SOURCE_FILES})