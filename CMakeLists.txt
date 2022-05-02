cmake_minimum_required(VERSION 3.16)

project(pro)

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

message(${OpenCV_INCLUDE_DIRS})

aux_source_directory(. DIR_SRCS)

SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

add_executable(pro ${DIR_SRCS})

target_link_libraries(pro ${OpenCV_LIBS})

