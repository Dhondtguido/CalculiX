include_guard(GLOBAL)

# user overrides
set(ARPACK_ROOT "" CACHE PATH "Prefix where ARPACK is installed")
set(ARPACK_LIBRARY_DIR "" CACHE FILEPATH "Full path to libarpack")

# resolve library
if(ARPACK_LIBRARY_DIR)

  if(NOT EXISTS "${ARPACK_LIBRARY_DIR}")
    message(FATAL_ERROR "ARPACK_LIBRARY_DIR does not exist: ${ARPACK_LIBRARY_DIR}")
  endif()

  set(_arpack_lib "${ARPACK_LIBRARY_DIR}")

else()

  find_library(_arpack_lib
    NAMES arpack arpack_ng
    HINTS
    ${ARPACK_ROOT}
    PATH_SUFFIXES
    lib lib64 lib/x86_64-linux-gnu
  )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(ARPACK
  REQUIRED_VARS _arpack_lib
)

if(NOT ARPACK_FOUND)
  message(FATAL_ERROR
    "ARPACK not found.\n"
    "Please provide either:\n"
    "  -DARPACK_LIBRARY_DIR=/path/to/libarpack\n"
    "or\n"
    "  -DARPACK_ROOT=/path/to/arpack/root\n"
  )
endif()

if(ARPACK_FOUND AND NOT TARGET ARPACK::ARPACK)

  add_library(ARPACK::ARPACK UNKNOWN IMPORTED)

  set_target_properties(ARPACK::ARPACK PROPERTIES
    IMPORTED_LOCATION "${_arpack_lib}"
  )

endif()
