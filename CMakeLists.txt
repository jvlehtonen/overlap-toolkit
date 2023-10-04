cmake_minimum_required(VERSION 3.10)

project(overlap LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)
set(CMAKE_BUILD_TYPE Release)

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "~/.local")
endif()

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Core REQUIRED)

add_executable(o-lap src/o-lap.cpp src/Mol2Read.cpp src/Point.cpp src/Atom.cpp src/json.cpp)

target_link_libraries(o-lap Qt5::Core)

install(TARGETS o-lap DESTINATION bin)
install(FILES data/cutoffs.json data/atomtypes.json DESTINATION share/SBL/o-lap)
