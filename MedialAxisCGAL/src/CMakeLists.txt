cmake_minimum_required(VERSION 3.12)
project(MedialAxisIpelet LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(CGAL REQUIRED)

get_filename_component(target_name MedialAxisIpelet.cpp NAME_WE)

add_executable(${target_name}
  MedialAxis.cpp
  MedialAxisIpelet.cpp
)

target_link_libraries(${target_name} PRIVATE CGAL::CGAL)
