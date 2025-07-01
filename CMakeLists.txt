cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

include(ExternalProject)

project(TexNeutSim)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Locate Geant4 package
find_package(Geant4 REQUIRED ui_all vis_all)

# Locate ROOT package
find_package(ROOT 6 QUIET CONFIG REQUIRED COMPONENTS RIO Tree Hist)

include(${Geant4_USE_FILE})

set(PROJECT_SOURCE_DIR_GEANT ${PROJECT_SOURCE_DIR}/src/)

file(GLOB sources ${PROJECT_SOURCE_DIR_GEANT}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR_GEANT}/*.hh)

file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

file(GLOB MACRO_FILES "*.dat")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

# Define dictionary file path
set(PCMFILE ${CMAKE_BINARY_DIR}/libLI6SIM_rdict.pcm)
add_definitions(-DPCMFILE=\"${PCMFILE}\")

# Define shared library file path (this also forces the correct file name)
set(SOFILE ${CMAKE_BINARY_DIR}/libLI6SIM.so)
add_definitions(-DSOFILE=\"${SOFILE}\")

# Set simlib directory
set(SIMLIB /home/Li6Webb/simlib)

# Include the MENATE_R library
add_subdirectory(lib/MENATE_R)

# Include the charged particle simulation library
ExternalProject_Add(li6sim
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/lib/Li6sim                          # Path to the source directory of the external project
	BINARY_DIR ${CMAKE_BINARY_DIR}/lib/Li6sim-build                    # Path to build directory
	CMAKE_ARGS
		-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/lib/Li6sim-install    # Installation path
		-DLI6SIM_LIBRARY_ONLY=ON                                         # Build the library instead of the executables
	BUILD_COMMAND ${CMAKE_COMMAND} --build . --target LI6SIM           # Custom build target
	INSTALL_COMMAND ""                                                 # Skip the installation step
	BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/lib/Li6sim-build/libLI6SIM.so # Path to built library
)

add_library(LI6SIM_IMPORTED STATIC IMPORTED GLOBAL)
set_target_properties(LI6SIM_IMPORTED PROPERTIES
	IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/Li6sim-build/libLI6SIM.so
	INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib/Li6sim/src;${SIMLIB}"
)

# Ensure the external project is built before this target is used
add_dependencies(LI6SIM_IMPORTED li6sim)

add_executable(texneutsim ${PROJECT_SOURCE_DIR_GEANT}/texneutsim.cc ${sources} ${headers})
target_link_libraries(texneutsim ${Geant4_LIBRARIES} MENATE_R LI6SIM_IMPORTED ROOT::RIO ROOT::Tree ROOT::Hist)
set_target_properties(texneutsim PROPERTIES
	BUILD_RPATH "${CMAKE_BINARY_DIR}"
	INSTALL_RPATH "${CMAKE_BINARY_DIR}"
)

# Set up copying of ROOT dictionary files and library after project is built
ExternalProject_Get_Property(li6sim BINARY_DIR)

add_custom_command(
	OUTPUT ${PCMFILE} ${CMAKE_BINARY_DIR}/libLI6SIM.rootmap ${SOFILE}
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		${BINARY_DIR}/libLI6SIM_rdict.pcm
		${PCMFILE}
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		${BINARY_DIR}/libLI6SIM.rootmap
		${CMAKE_BINARY_DIR}/libLI6SIM.rootmap
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		${BINARY_DIR}/libLI6SIM.so
		${SOFILE}
	DEPENDS li6sim
	COMMENT "Copying ROOT dictionary files and library from li6sim build directory"
)

add_custom_target(copy_li6sim_dict_files ALL
    DEPENDS ${PCMFILE} ${CMAKE_BINARY_DIR}/libLI6SIM.rootmap ${SOFILE}
)

