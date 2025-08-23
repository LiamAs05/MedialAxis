cmake_minimum_required(VERSION 3.12...3.29)
project(MedialAxis)
set(CMAKE_CXX_STANDARD 20)

find_package(CGAL REQUIRED)

set(IPE_ROOT "C:/Users/liamd/Desktop/ipe-7.2.29")
set(IPE_INCLUDE_DIR "${IPE_ROOT}/include")
set(IPE_LIB_DIR     "${IPE_ROOT}/lib")
set(IPE_DLL_DIR     "${IPE_ROOT}/bin")

include_directories(${IPE_INCLUDE_DIR})
link_directories(${IPE_LIB_DIR})

file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

foreach(cppfile ${cppfiles})
  get_filename_component(target_name ${cppfile} NAME_WE)
  add_executable(${target_name} ${cppfile})
  target_link_libraries(${target_name} PRIVATE CGAL::CGAL ipe ipeui ipecanvas)

  add_custom_command(TARGET ${target_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${IPE_DLL_DIR}/ipe.dll"
      "${IPE_DLL_DIR}/ipeui.dll"
      "${IPE_DLL_DIR}/ipecanvas.dll"
      $<TARGET_FILE_DIR:${target_name}>
  )
endforeach()
