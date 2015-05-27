cmake_minimum_required(VERSION 2.8)
 
PROJECT(MoCo)

FIND_PACKAGE(ITK REQUIRED)

INCLUDE(${ITK_USE_FILE})

SET( VTK_DIR /Users/matheusviana/Desktop/VTK-BIN)
set( VTK_DIR ${VTK_DIR} CACHE FILEPATH "" FORCE)

message(STATUS "VTK_DIR='${VTK_DIR}'")

if (ITKVtkGlue_LOADED)

	find_package(VTK REQUIRED)
	include(${VTK_USE_FILE})


else()

  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)

endif()

ADD_EXECUTABLE(MoCo MACOSX_BUNDLE MoCo.cxx)

SET( EXT_LIBRARIES ${EXT_LIBRARIES} vtkImagingCore vtkCommonMisc vtkFiltersSources vtkFiltersCore vtkImagingFourier vtkImagingGeneral vtkCommonDataModel vtkCommonCore vtkCommonTransforms vtkIOLegacy vtkCommonExecutionModel vtkIOCore vtkIOImage vtkRenderingFreeType vtkFiltersGeneral)
TARGET_LINK_LIBRARIES(MoCo  ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${EXT_LIBRARIES})

message(STATUS "VTK_LIB='${VTK_LIBRARIES}'")
