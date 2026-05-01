#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "JKQ::DDpackage" for configuration "Release"
set_property(TARGET JKQ::DDpackage APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(JKQ::DDpackage PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libdd_package.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS JKQ::DDpackage )
list(APPEND _IMPORT_CHECK_FILES_FOR_JKQ::DDpackage "${_IMPORT_PREFIX}/lib/libdd_package.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
