include_guard(GLOBAL)

# user overrides
set(SPOOLES_ROOT "" CACHE PATH "Prefix where SPOOLES is installed")
set(SPOOLES_LIBRARY_DIR "" CACHE FILEPATH "Full path to libspooles")
set(SPOOLES_INCLUDE_DIR "" CACHE PATH "Include directory for SPOOLES headers")

# resolve library
if(SPOOLES_LIBRARY_DIR)

  if(NOT EXISTS "${SPOOLES_LIBRARY_DIR}")
    message(FATAL_ERROR "SPOOLES_LIBRARY_DIR does not exist: ${SPOOLES_LIBRARY_DIR}")
  endif()

  set(_spooles_lib "${SPOOLES_LIBRARY_DIR}")

else()

  find_library(_spooles_lib
    NAMES spooles
    HINTS
    ${SPOOLES_ROOT}
    PATH_SUFFIXES
    lib lib64 lib/x86_64-linux-gnu
  )

endif()

# resolve headers
find_path(_spooles_include_dir
  NAMES SPOOLES.h
  PATH_SUFFIXES spooles
)

if(NOT _spooles_include_dir AND SPOOLES_ROOT)
  find_path(_spooles_include_dir
    NAMES SPOOLES.h
    HINTS "${SPOOLES_ROOT}/include/spooles"
    PATH_SUFFIXES "*/spooles"
  )
endif()

if(NOT _spooles_include_dir AND SPOOLES_INCLUDE_DIR)
  if(EXISTS "${SPOOLES_INCLUDE_DIR}/SPOOLES.h")
    set(_spooles_include_dir "${SPOOLES_INCLUDE_DIR}")
  endif()
endif()

if(NOT _spooles_include_dir)
  message(FATAL_ERROR
    "SPOOLES headers not found (looking for SPOOLES.h).\n"
    "Please provide:\n"
    "  -DSPOOLES_INCLUDE_DIR=/path/to/spooles/headers\n"
  )
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(SPOOLES
  REQUIRED_VARS _spooles_lib _spooles_include_dir
)

if(NOT SPOOLES_FOUND)
  message(FATAL_ERROR
    "SPOOLES not found.\n"
    "Please provide either:\n"
    "  -DSPOOLES_LIBRARY_DIR=/path/to/libspooles\n"
    "or\n"
    "  -DSPOOLES_ROOT=/path/to/spooles/root\n"
    "Headers will be resolved from SPOOLES_ROOT or\n"
    "  -DSPOOLES_INCLUDE_DIR=/path/to/spooles/headers\n"
  )
endif()

if(SPOOLES_FOUND AND NOT TARGET SPOOLES::SPOOLES)

  add_library(SPOOLES::SPOOLES UNKNOWN IMPORTED)

  set_target_properties(SPOOLES::SPOOLES PROPERTIES
    IMPORTED_LOCATION "${_spooles_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_spooles_include_dir}"
  )

endif()
