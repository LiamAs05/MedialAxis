cmake_minimum_required(VERSION 3.15)
project(MedialAxisIpelet LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -----------------------------
# Find CGAL
# -----------------------------
find_package(CGAL REQUIRED)

# -----------------------------
# Ipe configuration
# -----------------------------
# Adjust path to your Ipe build or install
set(IPE_ROOT "/usr/local/")   # or wherever Ipe is installed
set(IPE_INCLUDE_DIR "${IPE_ROOT}/include/ipe")
set(IPE_LIB_DIR     "${IPE_ROOT}/lib/ipe/7.3.1/")

include_directories(${IPE_INCLUDE_DIR})
link_directories(${IPE_LIB_DIR})

# -----------------------------
# Build ipelet as shared module
# -----------------------------
add_library(medialaxis MODULE
    MedialAxis.cpp
    MedialAxisIpelet.cpp
)

# CGAL
target_link_libraries(medialaxis PRIVATE CGAL::CGAL)

# Ipe (order matters!)
target_link_libraries(medialaxis PRIVATE libipe.so)

# -----------------------------
# Output name and extension
# -----------------------------
# On Linux, Ipe expects `.so` or `.ipelet`
set_target_properties(medialaxis PROPERTIES
    PREFIX ""                                # no "lib" prefix
)

# -----------------------------
# Install or copy
# -----------------------------
# Option 1: install into Ipe's ipelets directory
install(TARGETS medialaxis
    DESTINATION "${IPE_ROOT}/lib/ipe/ipelets"
)

# Option 2: or copy next to your binary dir for manual testing
add_custom_command(TARGET medialaxis POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:medialaxis>
        ${CMAKE_BINARY_DIR}/medialaxis.ipelet
)

