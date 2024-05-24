# Try to find hwloc libraries and headers.
# Once done this will define
#  HWLOC_FOUND - System has hwloc libraries and headers
#  HWLOC_INCLUDE_DIRS - The hwloc include directory
#  HWLOC_LIBRARIES - The libraries needed to use hwloc

find_path(
  HWLOC_PREFIX
  NAMES include/hwloc.h
)

if (NOT HWLOC_PREFIX AND NOT $ENV{HWLOC_BASE} STREQUAL "")
  set(HWLOC_PREFIX $ENV{HWLOC_BASE})
endif()

find_path(HWLOC_INCLUDE_DIR NAMES hwloc.h
  HINTS
  ${HWLOC_PREFIX}/include
)

find_library(HWLOC_LIBRARIES NAMES hwloc
  HINTS
  ${HWLOC_PREFIX}/lib
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  HWLOC DEFAULT_MSG
  HWLOC_LIBRARIES
  HWLOC_INCLUDE_DIRS
)

mark_as_advanced(
  HWLOC_LIBRARIES
  HWLOC_INCLUDE_DIRS
)
