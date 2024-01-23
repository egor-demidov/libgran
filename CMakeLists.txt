cmake_minimum_required(VERSION 3.26)
project(libgran)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3 -flto=auto")

# SYSTEM keyword suppressed warnings from included headers
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)
include_directories(include)

find_package(TBB REQUIRED)

add_executable(libgran main.cpp writer.cpp)

target_link_libraries(libgran PRIVATE TBB::tbb)
