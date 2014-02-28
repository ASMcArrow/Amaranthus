cmake_minimum_required(VERSION 2.8)
 
project(Ocular)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if(ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
endif()

#set(CMAKE_AUTOMOC ON)
#find_package(Qt5Widgets)

add_executable(ImageToVTKImageFilter main.cpp)
target_link_libraries(ImageToVTKImageFilter
${ItkVtkGlue}
${VTK_LIBRARIES}
${ITK_LIBRARIES})

# Use the Widgets module from Qt 5.
#qt5_use_modules(ImageToVTKImageFilter Widgets)

