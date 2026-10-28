# Fallback find module for PaStiX4CalculiX, which ships no CMake config package.
#
# PaStiX4CalculiX is discovered via pkg-config.

include_guard(GLOBAL)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PASTIX4CALCULIX pastix)

if(NOT PASTIX4CALCULIX_FOUND)
  message(FATAL_ERROR
    "PaStiX4CalculiX not found via pkg-config (module 'pastix').\n"
    "Ensure pastix.pc is accessible from PKG_CONFIG_PATH: Add the directory\n"
    "containing pastix.pc (e.g. <prefix>/lib/pkgconfig) to PKG_CONFIG_PATH\n"
    ", then reconfigure."
  )
endif()

if(NOT TARGET PASTIX::pastix)
  add_library(PASTIX::pastix UNKNOWN IMPORTED)
  set(_PASTIX_LIB_LOC "")
  foreach(_dir IN LISTS PASTIX4CALCULIX_LIBRARY_DIRS)
    if(EXISTS "${_dir}/libpastix.so")
      set(_PASTIX_LIB_LOC "${_dir}/libpastix.so")
      break()
    elseif(EXISTS "${_dir}/libpastix.a")
      set(_PASTIX_LIB_LOC "${_dir}/libpastix.a")
      break()
    endif()
  endforeach()
  if(NOT _PASTIX_LIB_LOC)
    message(FATAL_ERROR "libpastix not found in: ${PASTIX4CALCULIX_LIBRARY_DIRS}")
  endif()
  set_target_properties(PASTIX::pastix PROPERTIES
    IMPORTED_LOCATION "${_PASTIX_LIB_LOC}"
    INTERFACE_INCLUDE_DIRECTORIES "${PASTIX4CALCULIX_INCLUDE_DIRS}"
    INTERFACE_LINK_DIRECTORIES "${PASTIX4CALCULIX_LIBRARY_DIRS}"
    INTERFACE_LINK_LIBRARIES "${PASTIX4CALCULIX_LIBRARIES}"
  )
endif()
