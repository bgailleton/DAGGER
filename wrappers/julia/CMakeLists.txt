project(JuDagger)

cmake_minimum_required(VERSION 3.5)
set(CMAKE_MACOSX_RPATH 1)
add_compile_definitions(DAGGER_FT_JULIA)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

add_library(judagger SHARED main.cpp)

include_directories("../../DAGGER/")

target_link_libraries(judagger JlCxx::cxxwrap_julia)

install(TARGETS
  judagger
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib)
