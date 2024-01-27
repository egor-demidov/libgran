cmake_minimum_required(VERSION 3.22)
project(libgran)

set(CMAKE_CXX_STANDARD 17)

include(CTest)

set(ERROR_ON_WARN OFF)

if (${CMAKE_COMPILER_IS_GNUCXX})
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native -flto=auto")
    if (${ERROR_ON_WARN})
        set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")
    endif ()
endif ()

# SYSTEM keyword suppressed warnings from included headers
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)
include_directories(include)

if (${CMAKE_COMPILER_IS_GNUCXX})
    find_package(TBB REQUIRED)
endif ()

add_executable(libgran main.cpp writer.cpp)
# A particle and a lattice colliding
add_executable(contact_test test/contact.cpp test/compute_energy.cpp)
# A particle and a lattice with VdW interactions colliding
add_executable(hamaker_test test/hamaker.cpp writer.cpp)
# Two cubes with VdW interactions colliding
add_executable(hamaker_2_test test/hamaker_2.cpp writer.cpp)
# A conical lattice and an octahedral lattice colliding
add_executable(hamaker_3_test test/hamaker_3.cpp writer.cpp)
# A sintered lattice undergoing rigid-body-motion
add_executable(sintered_test test/sintered.cpp writer.cpp)

if (${CMAKE_COMPILER_IS_GNUCXX})
    target_link_libraries(libgran PRIVATE TBB::tbb)
    target_link_libraries(contact_test PRIVATE TBB::tbb)
    target_link_libraries(hamaker_test PRIVATE TBB::tbb)
    target_link_libraries(hamaker_2_test PRIVATE TBB::tbb)
    target_link_libraries(hamaker_3_test PRIVATE TBB::tbb)
    target_link_libraries(sintered_test PRIVATE TBB::tbb)
endif ()

add_test(NAME contact_test COMMAND ${CMAKE_BINARY_DIR}/contact_test)
