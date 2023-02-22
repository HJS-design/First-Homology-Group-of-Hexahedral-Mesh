# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
cmake_minimum_required(VERSION 3.1...3.15)
project( Surface_mesh_Examples )

#CGAL_Qt5 is needed for the drawing.
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(libigl  REQUIRED)


find_package(CGAL REQUIRED)

create_single_source_cgal_program("draw_surface_mesh.cpp")
target_link_libraries(draw_surface_mesh PUBLIC  igl::core igl::common igl::opengl igl::opengl_glfw  )
#end of the file
