cmake_minimum_required(VERSION 3.16)
project(MedialAxis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "Release")

# ---- Options ----
option(BUILD_CLI "Build standalone CLI tester" ON)
option(BUILD_IPELET "Build Ipelet plugin" ON)

# ---- CGAL ----
find_package(CGAL REQUIRED)

# ---- Ipe (only for Ipelet) ----
if(BUILD_IPELET)
  if(NOT DEFINED IPE_PREFIX)
    message(FATAL_ERROR "Please pass -DIPE_PREFIX=/path/to/ipe/install")
  endif()

  set(IPE_INCLUDE_DIR "${IPE_PREFIX}/include")
  set(IPE_LIB_DIR "${IPE_PREFIX}/lib")

  include_directories(${IPE_INCLUDE_DIR})
  link_directories(${IPE_LIB_DIR})
endif()

# ---- Sources ----
set(CORE_SRCS
  MedialAxis.cpp
  MedialAxis.hpp
)

# ---- CLI executable ----
if(BUILD_CLI)
  add_executable(medial_axis_cli
    cli_main.cpp
    ${CORE_SRCS}
  )
  target_link_libraries(medial_axis_cli PRIVATE CGAL::CGAL)
endif()

# ---- Ipelet plugin ----
if(BUILD_IPELET)
  add_library(medialaxis SHARED
    MedialAxisIpelet.cpp
    ${CORE_SRCS}
  )
  target_link_libraries(medialaxis PRIVATE CGAL::CGAL libipe.so)

  set_target_properties(medialaxis PROPERTIES
    PREFIX "" # removes "lib" prefix, so it builds medialaxis_ipelet.so not libmedialaxis_ipelet.so
  )
endif()

