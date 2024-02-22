cmake_minimum_required(VERSION 3.22)
project(libgran)

set(CMAKE_CXX_STANDARD 17)

include(CTest)

set(ERROR_ON_WARN OFF)

set(CMAKE_LIBGRAN_USE_OMP ON)

if (${CMAKE_LIBGRAN_USE_OMP})
    add_compile_definitions(LIBGRAN_USE_OMP)
endif ()

if (${APPLE})
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -fopenmp -flto -march=native")
    if (${ERROR_ON_WARN})
        set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")
    endif ()
    set(LIBRARY_LIST "-lomp -flto")
endif ()

if (${MSVC})
    add_compile_definitions(_USE_MATH_DEFINES)
    set(CMAKE_CXX_FLAGS "/O2 /openmp /EHsc /GL /fp:except /Wall")
    if (${ERROR_ON_WARN})
        set(CMAKE_CXX_FLAGS "/Wx ${CMAKE_CXX_FLAGS}")
    endif ()
    set(LIBRARY_LIST "")
endif ()

if (${CMAKE_COMPILER_IS_GNUCXX})
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native -flto=auto -fopenmp")
    if (${ERROR_ON_WARN})
        set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")
    endif ()
    find_package(TBB REQUIRED)
    set(LIBRARY_LIST PRIVATE TBB::tbb)
endif ()

# SYSTEM keyword suppressed warnings from included headers
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)
include_directories(include)

add_executable(libgran main.cpp writer.cpp test/compute_energy.cpp)
# A particle and a lattice colliding
add_executable(contact_test test/contact.cpp test/compute_energy.cpp)
# A particle and a lattice with VdW interactions colliding
add_executable(hamaker_test test/hamaker.cpp test/compute_energy.cpp test/mass_distribution.cpp)
# Two cubes with VdW interactions colliding
add_executable(hamaker_2_test test/hamaker_2.cpp writer.cpp)
# A conical lattice and an octahedral lattice colliding
add_executable(hamaker_3_test test/hamaker_3.cpp writer.cpp)
# A sintered lattice undergoing rigid-body-motion
add_executable(sintered_test test/sintered.cpp test/mass_distribution.cpp)
# A sintered lattice undergoing rigid-body motion with the new sintering model
add_executable(alt_sintered_test test/alt_sintered.cpp writer.cpp test/mass_distribution.cpp)
# Particle is collided with a plane
add_executable(surface_test test/surface.cpp writer.cpp test/compute_energy.cpp)

add_executable(sintered_2_test test/sintered_2.cpp writer.cpp)

target_link_libraries(libgran ${LIBRARY_LIST})
target_link_libraries(contact_test ${LIBRARY_LIST})
target_link_libraries(hamaker_test ${LIBRARY_LIST})
target_link_libraries(hamaker_2_test ${LIBRARY_LIST})
target_link_libraries(hamaker_3_test ${LIBRARY_LIST})
target_link_libraries(sintered_test ${LIBRARY_LIST})
target_link_libraries(sintered_2_test ${LIBRARY_LIST})
target_link_libraries(alt_sintered_test ${LIBRARY_LIST})
target_link_libraries(surface_test ${LIBRARY_LIST})

add_test(NAME contact_test COMMAND ${CMAKE_BINARY_DIR}/contact_test)
add_test(NAME hamaker_test COMMAND ${CMAKE_BINARY_DIR}/hamaker_test)
add_test(NAME sintered_test COMMAND ${CMAKE_BINARY_DIR}/sintered_test)
add_test(NAME alt_sintered_test COMMAND ${CMAKE_BINARY_DIR}/alt_sintered_test)
add_test(NAME surface_test COMMAND ${CMAKE_BINARY_DIR}/surface_test)
