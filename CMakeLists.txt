cmake_minimum_required(VERSION 3.26)
project(libgran)

set(CMAKE_CXX_STANDARD 17)

include(CTest)

set(ERROR_ON_WARN OFF)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native -flto=auto")
if (${ERROR_ON_WARN})
    set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")
endif ()

# SYSTEM keyword suppressed warnings from included headers
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)
include_directories(include)

find_package(TBB REQUIRED)

add_executable(libgran main.cpp writer.cpp)
add_executable(contact_test test/contact.cpp writer.cpp)

target_link_libraries(libgran PRIVATE TBB::tbb)
target_link_libraries(contact_test PRIVATE TBB::tbb)
