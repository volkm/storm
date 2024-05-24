# - Try to find libz3
# Once done this will define
#  LIBZ3_FOUND - System has libz3
#  LIBZ3_INCLUDE_DIRS - The libz3 include directories
#  LIBZ3_LIBRARIES - The libraries needed to use libz3

find_path(Z3_INCLUDE_DIRS NAMES z3++.h
    PATHS ENV PATH INCLUDE "${Z3_ROOT}/include" "/usr/include/z3" "/usr/local/include/z3/"
)

find_library(Z3_LIBRARIES NAMES z3
    PATHS ENV PATH INCLUDE "${Z3_ROOT}/lib"
)
if(Z3_INCLUDE_DIRS AND EXISTS "${Z3_INCLUDE_DIRS}/z3_version.h")
    file(STRINGS ${Z3_INCLUDE_DIRS}/z3_version.h Z3_VERSION_MAJOR REGEX "^#define[\t ]+Z3_MAJOR_VERSION .*")
    file(STRINGS ${Z3_INCLUDE_DIRS}/z3_version.h Z3_VERSION_MINOR REGEX "^#define[\t ]+Z3_MINOR_VERSION .*")
    file(STRINGS ${Z3_INCLUDE_DIRS}/z3_version.h Z3_VERSION_PATCH REGEX "^#define[\t ]+Z3_BUILD_NUMBER .*")
    string(REGEX MATCH "[0-9]+$" Z3_VERSION_MAJOR "${Z3_VERSION_MAJOR}")
    string(REGEX MATCH "[0-9]+$" Z3_VERSION_MINOR "${Z3_VERSION_MINOR}")
    string(REGEX MATCH "[0-9]+$" Z3_VERSION_PATCH "${Z3_VERSION_PATCH}")
    set(Z3_VERSION "${Z3_VERSION_MAJOR}.${Z3_VERSION_MINOR}.${Z3_VERSION_PATCH}")
endif()

# handle the QUIETLY and REQUIRED arguments and set SPOT_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Z3
	                          DEFAULT_MSG
				  REQUIRED_VARS Z3_LIBRARIES Z3_INCLUDE_DIRS
	                          VERSION_VAR Z3_VERSION
)


mark_as_advanced(Z3_LIBRARIES Z3_INCLUDE_DIRS)
