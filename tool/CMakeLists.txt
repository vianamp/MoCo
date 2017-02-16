
cmake_minimum_required(VERSION 2.8)

PROJECT(CreatePolyData)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(CreatePolyData MACOSX_BUNDLE CreatePolyData)

if(VTK_LIBRARIES)
  target_link_libraries(CreatePolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(CreatePolyData vtkHybrid vtkWidgets)
endif()
