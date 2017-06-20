cmake_minimum_required(VERSION 3.4)

enable_language(Fortran)
project(polysample)

# For a large number of source files you can create it in a simpler way
# using file() function:
set(SRC polysample.f90)

include_directories(${CMAKE_SOURCE_DIR}/../include)

# Create an library from sources
add_library(polysample ${SRC})

enable_testing()
add_subdirectory(test)

add_test(NAME test1 COMMAND test1)
