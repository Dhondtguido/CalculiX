# =============================================================================
# CMake project: CalculiX
#
# Changelog:
#   Date        Change                        Author
#   2026-10-08  Introduced CMake              Christoph WOELFLE
# =============================================================================

# patching hooks
# CCX_SOURCE_EXCLUDES       : list of src/<basename> to remove (e.g. "utemp.f")
# CCX_SOURCE_ADDITIONS      : absolute paths of extra source files to compile
# CCX_EXTRA_COMPILE_DEFS    : additional compile definitions
# CCX_EXTRA_COMPILE_OPTIONS : additional compile options
# CCX_EXTRA_INCLUDE         : additional include directories
# CCX_EXTRA_LIBS            : additional link libraries

cmake_minimum_required(VERSION 3.28)

project(CalculiX)
enable_language(C Fortran)

# options
option(CALCULIX_INT64 "build CalculiX with 64-bit integers (ILP64)" OFF)
option(SPOOLES "enable SPOOLES solver" ON)
option(PARDISO "enable Pardiso solver" OFF)
option(PASTIX "enable PaStiX solver" OFF)
option(PASTIX4CALCULIX "enable forked PaStiX4CalculiX solver" OFF)
option(MKL "use Intel OneAPI MKL (BLAS/LAPACK)" OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FLAGS)

# glob sources
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.f"
)
# remove entry point and includes
list(REMOVE_ITEM SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CalculiX.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/gauss.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/xlocal.f"
)
# PaStiX v6.4.0 compatibility: substitute src/pastix.c
if(PASTIX)
  list(APPEND CCX_SOURCE_EXCLUDES "pastix.c")
  list(APPEND CCX_SOURCE_ADDITIONS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compat/pastix/v6.4.0/pastix.c")
endif()

# process source patching hooks
foreach(_f IN LISTS CCX_SOURCE_EXCLUDES)
  list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/${_f}")
endforeach()
list(APPEND SOURCES ${CCX_SOURCE_ADDITIONS})

# machinery for updating of build time info in the binaries
set(CCX_BUILD_TIMESTAMP "${CMAKE_CURRENT_BINARY_DIR}/build.timestamp")
set(CCX_BUILD_TIMESTAMP_DEPS "${CMAKE_CURRENT_SOURCE_DIR}/src/CalculiX.c" ${SOURCES})
list(REMOVE_ITEM CCX_BUILD_TIMESTAMP_DEPS "${CMAKE_CURRENT_SOURCE_DIR}/src/buildinfo.c")
add_custom_command(
  OUTPUT "${CCX_BUILD_TIMESTAMP}"
  COMMAND ${CMAKE_COMMAND} -E touch ${CCX_BUILD_TIMESTAMP}
  DEPENDS ${CCX_BUILD_TIMESTAMP_DEPS}
)
add_custom_target(ccx_buildtime ALL DEPENDS "${CCX_BUILD_TIMESTAMP}")
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/buildinfo.c"
  PROPERTIES OBJECT_DEPENDS "${CCX_BUILD_TIMESTAMP}"
)

add_library(ccx ${SOURCES})
target_compile_definitions(ccx PUBLIC
  ARCH=${CMAKE_SYSTEM_NAME}
  MATRIXSTORAGE
)
target_link_libraries(ccx PUBLIC flags)

add_executable(CalculiX src/CalculiX.c)
target_link_libraries(CalculiX PRIVATE ccx ${CCX_EXTRA_LIBS} ccx)  # ccx twice on link line so CCX_EXTRA_LIB can use any symbol from ccx

# process compilation and linking patching hooks
if(CCX_EXTRA_COMPILE_DEFS)
  target_compile_definitions(ccx PUBLIC ${CCX_EXTRA_COMPILE_DEFS})
endif()
if(CCX_EXTRA_INCLUDE)
  target_include_directories(ccx PUBLIC ${CCX_EXTRA_INCLUDE})
endif()
if(CCX_EXTRA_COMPILE_OPTIONS)
  target_compile_options(ccx PUBLIC ${CCX_EXTRA_COMPILE_OPTIONS})
endif()

# integer type
if(CALCULIX_INT64)
  message(STATUS "CalculiX integer size: 64-bit (ILP64)")
  target_compile_definitions(ccx PUBLIC INTSIZE64)
else()
  message(STATUS "CalculiX integer size: 32-bit (LP64)")
endif()

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(ccx PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)

# BLAS/LAPACK
if(MKL)
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    message(FATAL_ERROR "MKL requires an Intel oneAPI compiler: its OpenMP runtime (libiomp5) "
      "must match the compiler's, otherwise two OpenMP runtimes get linked into CalculiX.")
  endif()
  if(CALCULIX_INT64)
    set(MKL_INTERFACE ilp64)
  else()
    set(MKL_INTERFACE lp64)
  endif()
  find_package(MKL CONFIG REQUIRED)
  message(STATUS "Found Intel MKL: ${MKL_ROOT}")
  target_link_libraries(ccx PUBLIC MKL::MKL)
else()
  if(CALCULIX_INT64)
    message(WARNING "building with CALCULIX_INT64=ON but MKL not found\n assert BLAS/LAPACK is ILP64 to avoid ABI mismatch")
  endif()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(ccx PUBLIC
    BLAS::BLAS
    LAPACK::LAPACK
  )
endif()

# ARPACK
find_package(arpackng CONFIG)
if(NOT arpackng_FOUND)
  find_package(ARPACK CONFIG)
  if(NOT ARPACK_FOUND)
    find_package(ARPACK MODULE REQUIRED)
  endif()
endif()
target_link_libraries(ccx PUBLIC ARPACK::ARPACK)
target_compile_definitions(ccx PUBLIC ARPACK)

# SPOOLES
if(SPOOLES)
  find_package(SPOOLES CONFIG QUIET)
  if(NOT SPOOLES_FOUND)
    find_package(SPOOLES MODULE REQUIRED)
  endif()
  target_link_libraries(ccx PUBLIC SPOOLES::SPOOLES)
  target_compile_definitions(ccx PUBLIC SPOOLES)
endif()

# PARDISO
if(PARDISO)
  if(NOT MKL_FOUND)
    message(FATAL_ERROR "MKL is required for builds with PARDISO")
  else()
    target_compile_definitions(ccx PRIVATE PARDISO)
  endif()
endif()

# PaStiX
if(PASTIX OR PASTIX4CALCULIX)
  if(PASTIX AND PASTIX4CALCULIX)
    message(FATAL_ERROR
      "PaStiX and PaStiX4CalculiX are mutually exclusive."
    )
  endif()
  if(NOT CALCULIX_INT64)
    message(WARNING
      "PaStiX with CALCULIX_INT64=OFF is an untested combination.")
  endif()
  if(PASTIX)
    find_package(PASTIX CONFIG REQUIRED)
  elseif(PASTIX4CALCULIX)
    find_package(PASTIX4CALCULIX MODULE REQUIRED)
  endif()
  target_link_libraries(ccx PUBLIC PASTIX::pastix)
  target_compile_definitions(ccx PUBLIC PASTIX)
endif()
