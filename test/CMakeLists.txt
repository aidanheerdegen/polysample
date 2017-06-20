cmake_minimum_required(VERSION 3.4)

enable_language(Fortran)

set(SRC test1.f90)

include_directories(${CMAKE_BINARY_DIR})

# Create an executable
add_executable(test1 ${SRC})

target_link_libraries(test1 polysample)

add_test(NAME test1 COMMAND test1)
