
# CMake Stuff
cmake_minimum_required(VERSION 3.20)

# must be set before project(...) call; version module is needed before
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(pythiaBeamEffects VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT (required for histograms and such)
find_package(ROOT REQUIRED COMPONENTS Core Hist RIO MathCore)

# Find HepMC3 (required because this is what the program is for)
find_package(HepMC3 REQUIRED)

# Find Pythia8
find_package(Pythia8 REQUIRED)
set(PYTHIA8_INCLUDE_DIRS ${Pythia8_INCLUDE_DIR})
set(PYTHIA8_LIBRARIES ${Pythia8_LIBRARY})
message(STATUS "Pythia8 found via FindPythia8.cmake")
message(STATUS "Pythia8 Include: ${PYTHIA8_INCLUDE_DIRS}")
message(STATUS "Pythia8 Libs: ${PYTHIA8_LIBRARIES}")

# Find Fastjet
find_package(fastjet REQUIRED)
set(FASTJET_INCLUDE_DIRS ${fastjet_INCLUDE_DIR})
set(FASTJET_LIBRARIES ${fastjet_LIBRARY})
message(STATUS "FastJet found via Findfastjet.cmake")
message(STATUS "FastJet Include: ${FASTJET_INCLUDE_DIRS}")
message(STATUS "FastJet Libs: ${FASTJET_LIBRARIES}")

# Set Output Directories
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib )

# Set Executable Name
set( PYTHIA_EXE runBeamShape )
set( PYTHIAHEPMC_EXE runBeamShapeHepMC )
set( BEAM_LIB eicBeamShape )

# Set Sources
set( PYTHIA_SRC PythiaBeamShape.cxx )
set( PYTHIAHEPMC_SRC PythiaBeamShapeHepMC.cxx )
set( BEAM_SRC eicBeamShape.cxx )

# Create library and executables
add_library( ${BEAM_LIB} ${BEAM_SRC} )
target_include_directories( ${BEAM_LIB} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHIA8_INCLUDE_DIRS} )
target_link_libraries( ${BEAM_LIB} PRIVATE ${PYTHIA8_LIBRARIES} ROOT::Core ROOT::Hist ROOT::RIO ROOT::MathCore ${ROOT_LIBRARIES} )

add_executable( ${PYTHIA_EXE} ${PYTHIA_SRC} )
target_include_directories( ${PYTHIA_EXE} PRIVATE ${FASTJET_INCLUDE_DIRS} )
target_link_libraries( ${PYTHIA_EXE} PRIVATE ${FASTJET_LIBRARIES} HepMC3::HepMC3 HepMC3::rootIO ROOT::Core ROOT::Hist ROOT::RIO ROOT::MathCore ${ROOT_LIBRARIES} ${BEAM_LIB} )

add_executable( ${PYTHIAHEPMC_EXE} ${PYTHIAHEPMC_SRC} )
target_link_libraries( ${PYTHIAHEPMC_EXE} PRIVATE HepMC3::HepMC3 HepMC3::rootIO ROOT::Core ROOT::Hist ROOT::RIO ROOT::MathCore ${ROOT_LIBRARIES} ${BEAM_LIB} )

install(TARGETS ${BEAM_LIB} LIBRARY)
install(TARGETS ${PYTHIA_EXE} ${PYTHIAHEPMC_EXE} RUNTIME)
