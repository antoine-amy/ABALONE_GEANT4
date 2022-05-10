cmake_minimum_required(VERSION 3.0.0)
project(ABALONE_GEANT4 VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(ABALONE_GEANT4 main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
