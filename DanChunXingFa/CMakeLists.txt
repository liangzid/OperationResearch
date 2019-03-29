cmake_minimum_required(VERSION 3.10.2)
project(DanChunXingFa)

set(CMAKE_CXX_STANDARD 14)

INCLUDE_DIRECTORIES(include)

AUX_SOURCE_DIRECTORY(src SRC_FILES)

add_executable(DanChunXingFa demo.cpp src include )