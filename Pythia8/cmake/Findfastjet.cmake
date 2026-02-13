# Find the FastJet includes and libraries.
#
# This module defines the `FastJet` imported target that encodes all
# necessary information in its target properties.

include(FindPackageHandleStandardArgs)

find_library(fastjet_LIBRARY
  NAMES fastjet
  PATHS /opt/local/lib 
  DOC "The fastjet library")
find_path(fastjet_INCLUDE_DIR
  NAMES fastjet
  PATHS /opt/local/include
  DOC "The fastjet include directory")

find_package_handle_standard_args(fastjet
  REQUIRED_VARS fastjet_LIBRARY fastjet_INCLUDE_DIR)

add_library(fastjet SHARED IMPORTED)
set_property(TARGET fastjet PROPERTY IMPORTED_LOCATION ${fastjet_LIBRARY})
set_property(TARGET fastjet PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${fastjet_INCLUDE_DIR})

mark_as_advanced(fastjet_FOUND fastjet_INCLUDE_DIR fastjet_LIBRARY)
