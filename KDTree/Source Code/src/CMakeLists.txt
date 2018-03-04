cmake_minimum_required(VERSION 2.8.11)

project(GeoLab-2016)

### Warning options #####################################################################
option(DISABLE_WARNINGS "Disable Warnings MSVC/GCC" ON)
if(DISABLE_WARNINGS)
 set(MSVC_WARN_LEVEL /W0)
 set(GCC_WARN_LEVEL -w)
else()
 set(MSVC_WARN_LEVEL /W3)
 set(GCC_WARN_LEVEL -Wall -Wno-long-long -pedantic)
endif(DISABLE_WARNINGS)
if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "${MSVC_WARN_LEVEL}" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSVC_WARN_LEVEL}")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_WARN_LEVEL} -std=c++11")
endif()
#########################################################################################

### Find VVRFramework ###################################################################
set(VVRFRAMEWORK_DIR CACHE PATH "Installation dir of VVRFramework")
link_directories(${VVRFRAMEWORK_DIR}/lib)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
#########################################################################################

### Specify Executables #################################################################
file(GLOB APP_SRC_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach(srcFile ${APP_SRC_FILES})
    string(REPLACE ".cpp" "" appName ${srcFile})
    if(EXISTS ${appName}.h)
        add_executable( ${appName} MACOSX_BUNDLE ${appName}.cpp ${appName}.h)
    else()
        add_executable( ${appName} MACOSX_BUNDLE ${appName}.cpp)
    endif()
    target_link_libraries(${appName} VVRScene GeoLib MathGeoLib ${Boost_LIBRARIES})
endforeach(srcFile ${APP_SRC_FILES})
#########################################################################################

### Set include path ####################################################################
include_directories(${VVRFRAMEWORK_DIR}/include)
include_directories(${VVRFRAMEWORK_DIR}/include/GeoLib)
include_directories(${VVRFRAMEWORK_DIR}/include/MathGeoLib)
#########################################################################################
