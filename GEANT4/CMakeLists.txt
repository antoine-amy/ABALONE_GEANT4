cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(AB)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

option(WITH_GEANT4_UIVIS "Build simulation with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(externals/g4sipm)
include_directories(externals/g4sipm/include)
include_directories(externals/g4sipm/externals/jansson/src)
include_directories(externals/g4sipm/externals/gtest/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Locate rnd sources and headers
file(GLOB sources_rnd ${PROJECT_SOURCE_DIR}/src/rnd/*.cc)
set(sources ${sources} ${sources_rnd})
file(GLOB headers_rnd ${PROJECT_SOURCE_DIR}/include/rnd/*.hh)
set(headers ${headers} ${headers_rnd})

# Add persistency sources and headers.
set(sources ${sources} "${PROJECT_SOURCE_DIR}/src/persistency/PersistencyHandler.cc" 
    "${PROJECT_SOURCE_DIR}/src/persistency/PersistencyHandlerMessenger.cc" 
    "${PROJECT_SOURCE_DIR}/src/persistency/PersistVisitor.cc" 
    "${PROJECT_SOURCE_DIR}/src/persistency/PersistVisitorFactory.cc" 
    "${PROJECT_SOURCE_DIR}/src/persistency/JsonPersistVisitor.cc")
set(headers ${headers} "${PROJECT_SOURCE_DIR}/include/persistency/PersistencyHandler.hh" 
    "${PROJECT_SOURCE_DIR}/include/persistency/PersistencyHandlerMessenger.hh" 
    "${PROJECT_SOURCE_DIR}/include/persistency/PersistVisitor.hh" 
    "${PROJECT_SOURCE_DIR}/include/persistency/PersistVisitorFactory.hh" 
    "${PROJECT_SOURCE_DIR}/include/persistency/JsonPersistVisitor.hh")
    
if(WITH_ROOT)
    set(sources ${sources} "${PROJECT_SOURCE_DIR}/src/persistency/RootPersistVisitor.cc")
    set(headers ${headers} "${PROJECT_SOURCE_DIR}/include/persistency/RootPersistVisitor.hh")
	set(CMAKE_CXX_FLAGS "-DWITH_ROOT ${CMAKE_CXX_FLAGS}")
endif()    

# Boost
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost 1.50.0 COMPONENTS date_time program_options filesystem system regex REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_subdirectory(plots)
add_subdirectory(resources)
add_subdirectory(run)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ABsimulation ABsimulation.cc ${sources} ${headers})
target_link_libraries(ABsimulation g4sipm boost_date_time jansson ${Geant4_LIBRARIES} ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build AB. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(ABSIMULATION_SCRIPTS
  ABsimulation.in
  ABsimulation.out
  full_run.mac
  beamon.mac
  execute.mac
  pe.mac
  decay.mac
  vis.mac
  gui.mac
  Biasdecay_z71.data
  )

foreach(_script ${ABSIMULATION_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# simulation standalone
#
add_custom_target(AB DEPENDS ABsimulation)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ABsimulation DESTINATION bin)
