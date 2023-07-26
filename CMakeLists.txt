cmake_minimum_required(VERSION 3.26)
project(loopSubdiv)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_USE_STATIC_LIBRARY       "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL              "Use OpenGL"         ON)
option(LIBIGL_WITH_PNG                 "Use PNG"            ON)
option(LIBIGL_WITH_OPENGL_GLFW         "Use GLFW"           ON)
find_package(LIBIGL REQUIRED)
# loop_subdiv
include_directories("include/")
file(GLOB SRCFILES src/*.cpp)

add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::opengl igl::opengl_glfw)

add_executable(subdivision "subdivision.cpp")
target_link_libraries(subdivision core igl::core igl::opengl igl::opengl_glfw)