add_custom_target(resources)
add_custom_target(test-resources)

set(STORM_3RDPARTY_SOURCE_DIR ${PROJECT_SOURCE_DIR}/resources/3rdparty)
set(STORM_3RDPARTY_BINARY_DIR ${PROJECT_BINARY_DIR}/resources/3rdparty)
set(STORM_3RDPARTY_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include/resources/3rdparty)

message(STATUS "Storm - Building external resources with ${STORM_RESOURCES_BUILD_JOBCOUNT} job(s) in parallel.")

#############################################################
##
##	l3pp
##
#############################################################

# Do not take a branch, needs internet connection.
ExternalProject_Add(
        l3pp_ext
        GIT_REPOSITORY https://github.com/hbruintjes/l3pp.git
        GIT_TAG e4f8d7fe6c328849aff34d2dfd6fd592c14070d5
        SOURCE_DIR ${STORM_3RDPARTY_SOURCE_DIR}/l3pp
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_INSTALL ON
)
ExternalProject_Get_Property(l3pp_ext source_dir)
set(l3pp_INCLUDE "${source_dir}/")
add_imported_library_interface(l3pp "${l3pp_INCLUDE}")
list(APPEND STORM_DEP_TARGETS l3pp)
add_dependencies(l3pp l3pp_ext)

#############################################################
##
##	gmm
##
#############################################################

add_imported_library_interface(gmm "${STORM_3RDPARTY_SOURCE_DIR}/gmm-5.2/include")
list(APPEND STORM_DEP_TARGETS gmm)

#############################################################
##
##	Eigen
##
#############################################################

# Checkout Eigen 
# We're checking out Eigen using the (currently not officially released) version 3.4.1
# See here https://gitlab.com/libeigen/eigen/-/commit/bae907b8f6078b1df290729eef946360315bd312
# The current official release is 3.4.0 which results in issues, see https://github.com/moves-rwth/storm/issues/162
message (STATUS "Storm - Including Eigen commit bae907b8f6078b1df290729eef946360315bd312 (unofficial 3.4.1 release).")
ExternalProject_Add(
        eigen_src
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG bae907b8f6078b1df290729eef946360315bd312
        SOURCE_DIR ${STORM_3RDPARTY_INCLUDE_DIR}/StormEigen
        PREFIX ${STORM_3RDPARTY_BINARY_DIR}/StormEigen-3.4.1alpha
        # First check whether patch was already applied (--reverse --check), otherwise apply patch
        PATCH_COMMAND git apply ${STORM_3RDPARTY_SOURCE_DIR}/patches/eigen341alpha.patch --reverse --check || git apply ${STORM_3RDPARTY_SOURCE_DIR}/patches/eigen341alpha.patch
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_INSTALL ON
)
add_imported_library_interface(StormEigen "${STORM_3RDPARTY_INCLUDE_DIR}")
list(APPEND STORM_DEP_TARGETS StormEigen)
add_dependencies(StormEigen eigen_src)



#############################################################
##
##	Boost
##
#############################################################

# Boost Option variables
set(Boost_USE_STATIC_LIBS ${USE_BOOST_STATIC_LIBRARIES})
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_NO_BOOST_CMAKE ON)

find_package(Boost 1.65.1 QUIET REQUIRED COMPONENTS filesystem system)
if (NOT Boost_FOUND)
    if (Boost_VERSION)
	message(FATAL_ERROR "The required Boost version is 1.65.1 or newer, however, only ${Boost_VERSION} was found.")
    else ()
   	message(FATAL_ERROR "Boost was not found.")
    endif ()
endif ()
if ((NOT Boost_LIBRARY_DIRS) OR ("${Boost_LIBRARY_DIRS}" STREQUAL ""))
    set(Boost_LIBRARY_DIRS "${Boost_INCLUDE_DIRS}/stage/lib")
endif ()

if (${Boost_VERSION} VERSION_GREATER_EQUAL "1.81.0")
    message(STATUS "Storm - Using workaround for Boost >= 1.81")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_PHOENIX_STL_TUPLE_H_")
endif()

set(CNTVAR 1)
foreach(BOOSTLIB ${Boost_LIBRARIES})
    add_imported_library(target-boost-${CNTVAR} SHARED ${BOOSTLIB} ${Boost_INCLUDE_DIRS})
    list(APPEND STORM_DEP_TARGETS target-boost-${CNTVAR}_SHARED)
    MATH(EXPR CNTVAR "${CNTVAR}+1")
endforeach()
message(STATUS "Storm - Using boost ${Boost_VERSION} (library version ${Boost_LIB_VERSION}).")
# set the information for the config header
set(STORM_BOOST_INCLUDE_DIR "${Boost_INCLUDE_DIRS}")

#############################################################
##
##	ExprTk
##
#############################################################

# Use the shipped version of ExprTK
message (STATUS "Storm - Including ExprTk.")
add_imported_library_interface(ExprTk "${STORM_3RDPARTY_SOURCE_DIR}/exprtk")
list(APPEND STORM_DEP_TARGETS ExprTk)



#############################################################
##
##	Parallel Hashmap
##
#############################################################

# Use the shipped version of Parallel Hashmap
message (STATUS "Storm - Including Parallel Hashmap.")
set(PHMAP_INCLUDE_DIR "${STORM_3RDPARTY_SOURCE_DIR}/parallel_hashmap/parallel_hashmap")
file(GLOB PHMAP_HEADERS "${PHMAP_INCLUDE_DIR}/*.h")

# Add the parallel_hashmap headers to the headers that are copied to the include directory in the build directory.
set(PHMAP_BINDIR_DIR ${STORM_3RDPARTY_INCLUDE_DIR}/parallel_hashmap)
include_directories("${PHMAP_BINDIR_DIR}")
foreach(HEADER ${PHMAP_HEADERS})
    string(REGEX REPLACE "${PROJECT_SOURCE_DIR}/src/?" "" RELATIVE_HEADER_PATH ${HEADER})
    string(REGEX MATCH "(.*)[/\\]" RELATIVE_DIRECTORY ${RELATIVE_HEADER_PATH})
    string(REGEX REPLACE "${RELATIVE_DIRECTORY}/?" "" HEADER_FILENAME ${RELATIVE_HEADER_PATH})
    add_custom_command(
        OUTPUT ${PHMAP_BINDIR_DIR}/parallel_hashmap/${HEADER_FILENAME}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${PHMAP_BINDIR_DIR}/parallel_hashmap
        COMMAND ${CMAKE_COMMAND} -E copy ${HEADER} ${PHMAP_BINDIR_DIR}/parallel_hashmap/${HEADER_FILENAME}
        DEPENDS ${PHMAP_INCLUDE_DIR}/${HEADER_FILENAME}
    )
    list(APPEND PHMAP_BINDIR_HEADERS ${PHMAP_BINDIR_DIR}/parallel_hashmap/${HEADER_FILENAME})
endforeach()

#############################################################
##
##	cpphoafparser
##
#############################################################

# Use the shipped version of cpphoafparser
message (STATUS "Storm - Including cpphoafparser 0.99.2")
include_directories("${STORM_3RDPARTY_SOURCE_DIR}/cpphoafparser-0.99.2/include")

#############################################################
##
##	ModernJSON
##
#############################################################

#use the shipped version of modernjson
message (STATUS "Storm - Including ModernJSON.")
add_imported_library_interface(ModernJSON "${STORM_3RDPARTY_SOURCE_DIR}/modernjson/include/")
list(APPEND STORM_DEP_TARGETS ModernJSON)

#############################################################
##
##	Z3 (optional)
##
#############################################################

find_package(Z3 QUIET)

# Z3 Defines
set(STORM_HAVE_Z3 ${Z3_FOUND})

if(Z3_FOUND)
    if(NOT EXISTS "${Z3_INCLUDE_DIR}/z3_version.h")
        message(FATAL_ERROR "No file z3_version.h found in ${Z3_INCLUDE_DIR}.")
    endif()

    # Parse z3 version from version file
    file(STRINGS ${Z3_INCLUDE_DIR}/z3_version.h Z3_VERSION_MAJOR REGEX "^#define[\t ]+Z3_MAJOR_VERSION .*")
    file(STRINGS ${Z3_INCLUDE_DIR}/z3_version.h Z3_VERSION_MINOR REGEX "^#define[\t ]+Z3_MINOR_VERSION .*")
    file(STRINGS ${Z3_INCLUDE_DIR}/z3_version.h Z3_VERSION_PATCH REGEX "^#define[\t ]+Z3_BUILD_NUMBER .*")
    string(REGEX MATCH "[0-9]+$" Z3_VERSION_MAJOR "${Z3_VERSION_MAJOR}")
    string(REGEX MATCH "[0-9]+$" Z3_VERSION_MINOR "${Z3_VERSION_MINOR}")
    string(REGEX MATCH "[0-9]+$" Z3_VERSION_PATCH "${Z3_VERSION_PATCH}")
    set(Z3_VERSION "${Z3_VERSION_MAJOR}.${Z3_VERSION_MINOR}.${Z3_VERSION_PATCH}")

    if(Z3_VERSION)
        # Split Z3 version into its components
        string(REPLACE "." ";" Z3_VERSION_LIST ${Z3_VERSION})
        list(GET Z3_VERSION_LIST 0 STORM_Z3_VERSION_MAJOR)
        list(GET Z3_VERSION_LIST 1 STORM_Z3_VERSION_MINOR)
        list(GET Z3_VERSION_LIST 2 STORM_Z3_VERSION_PATCH)

        if ("${Z3_VERSION}" VERSION_EQUAL "4.10.1")
            message (WARNING "Storm cannot be linked against Z3 version 4.10.1, see issue 252.")
            message(WARNING "Storm - Could not link with z3. Building of Prism/JANI models will not be supported.")
            set(Z3_FOUND FALSE)
        else()
            # Check whether the version of z3 supports optimization
            if(NOT "${Z3_VERSION}" VERSION_LESS "4.5.0")
                set(STORM_HAVE_Z3_OPTIMIZE ON)
                message (STATUS "Storm - Linking with Z3 (version ${Z3_VERSION}). (Z3 version supports optimization)")
            else()
                message (STATUS "Storm - Linking with Z3 (version ${Z3_VERSION}). (Z3 version does not support optimization)")
            endif()
            if (NOT "${Z3_VERSION}" VERSION_LESS "4.7.1")
                set(STORM_Z3_API_USES_STANDARD_INTEGERS ON)
            endif()

            add_imported_library(Z3 SHARED ${Z3_LIBRARIES} ${Z3_INCLUDE_DIRS})
            list(APPEND STORM_DEP_TARGETS Z3_SHARED)
        endif()
    else()
        message(WARNING "Storm - Could not obtain Z3 version. Building of Prism/JANI models will not be supported.")
        set(Z3_FOUND FALSE)
    endif()

else()
    message (WARNING "Storm - Z3 not found. Building of Prism/JANI models will not be supported.")
endif(Z3_FOUND)


#############################################################
##
##	glpk
##
#############################################################

include(${STORM_3RDPARTY_SOURCE_DIR}/include_glpk.cmake)

#############################################################
##
##	Gurobi (optional)
##
#############################################################

if (STORM_USE_GUROBI)
    find_package(GUROBI QUIET REQUIRED)
    set(STORM_HAVE_GUROBI ${GUROBI_FOUND})
    if (GUROBI_FOUND)
		if (EXISTS ${GUROBI_LIBRARY})
            message (STATUS "Storm - Linking with Gurobi (${GUROBI_CXX_LIBRARY}).")
            add_imported_library(GUROBI SHARED ${GUROBI_LIBRARY} ${GUROBI_INCLUDE_DIRS})
            list(APPEND STORM_DEP_TARGETS GUROBI_SHARED)
        else()
            # The FindGurobi.cmake script needs to be updated every now and then as the library file contains the version number...
            message(FATAL_ERROR "Gurobi Library ${GUROBI_LIBRARY} not found. If your Gurobi Version is higher then 9.0.0, please contact the Storm developers.")
        endif()
    else()
    	# The FindGurobi.cmake script needs to be updated every now and then as the library file contains the version number ...
        message(FATAL_ERROR "Gurobi Library requested but was not found. Make sure that GUROBI_ROOT points to the correct directory (containing include/ and lib/ subdirectories). If your Gurobi Version is higher then 9.0.0, please contact the Storm developers.")
    endif()
else()
    set(STORM_HAVE_GUROBI OFF)
endif()

#############################################################
##
##	CUDD
##
#############################################################
include(${STORM_3RDPARTY_SOURCE_DIR}/include_cudd.cmake)


#############################################################
##
##      carl
##
#############################################################

set(STORM_HAVE_CARL OFF)
set(CARL_MINVERSION "14.26")
set(CARL_C14VERSION "14")
if (NOT STORM_FORCE_SHIPPED_CARL)
    if (NOT "${STORM_CARL_DIR_HINT}" STREQUAL "")
		find_package(carl QUIET PATHS ${STORM_CARL_DIR_HINT} NO_DEFAULT_PATH)
    endif()
	if (NOT carl_FOUND)
		find_package(carl QUIET)
    endif()
endif()

set(STORM_SHIPPED_CARL OFF)
if(carl_FOUND AND NOT STORM_FORCE_SHIPPED_CARL)
    get_target_property(carlLOCATION lib_carl LOCATION)
    if("${carlLOCATION}" STREQUAL "carlLOCATION-NOTFOUND")
        if(EXISTS ${STORM_3RDPARTY_BINARY_DIR}/carl)
            message(WARNING "Storm - Library for carl location is not found but the directory ${STORM_3RDPARTY_BINARY_DIR}/carl exists. Will (re-)try to build a shipped version of carl.")
            set(STORM_SHIPPED_CARL ON)
        else()
            message(FATAL_ERROR "Library location for carl is not found, did you build carl?")
        endif()
    elseif(EXISTS ${carlLOCATION})
        #empty on purpose
    else()
        if(EXISTS ${STORM_3RDPARTY_BINARY_DIR}/carl)
            message(WARNING "Storm - File ${carlLOCATION} does not exist but the directory ${STORM_3RDPARTY_BINARY_DIR}/carl exists. Will (re-)try to build a shipped version of carl.")
            set(STORM_SHIPPED_CARL ON)
        else()
            message(FATAL_ERROR "File ${carlLOCATION} does not exist, did you build carl?")
        endif()
    endif()
    if("${carl_VERSION_MAJOR}" STREQUAL "${CARL_C14VERSION}")
        message(STATUS "Storm - Found carl-storm version")
        # empty on purpose. Maybe put a warning here?
        if("${carl_VERSION_MAJOR}.${carl_VERSION_MINOR}" VERSION_LESS "${CARL_MINVERSION}")
            message(FATAL_ERROR "Carl version outdated. We require ${CARL_MINVERSION}. Found ${carl_VERSION_MAJOR}.${carl_VERSION_MINOR} at ${carlLOCATION}")
        endif()
    else()
        message(FATAL_ERROR "We only support a diverged version of carl, indicated by Carl version 14.x. These versions can be found at https://github.com/moves-rwth/carl-storm. 
			      On this system, we found ${carl_VERSION_MAJOR}.${carl_VERSION_MINOR} at ${carlLOCATION}")
    endif()

    set(STORM_HAVE_CARL ON)
    message(STATUS "Storm - Use system version of carl.")
    message(STATUS "Storm - Linking with preinstalled carl ${carl_VERSION} (include: ${carl_INCLUDE_DIR}, library ${carl_LIBRARIES}, CARL_USE_CLN_NUMBERS: ${CARL_USE_CLN_NUMBERS}, CARL_USE_GINAC: ${CARL_USE_GINAC}).")
    set(STORM_HAVE_CLN ${CARL_USE_CLN_NUMBERS})
    set(STORM_HAVE_GINAC ${CARL_USE_GINAC})
else()
    set(STORM_SHIPPED_CARL ON)
endif()

if (STORM_SHIPPED_CARL)
	# The first external project will be built at *configure stage*
	message(STATUS "Carl - Start of config process")
	file(MAKE_DIRECTORY ${STORM_3RDPARTY_BINARY_DIR}/carl_download)
	execute_process(
	COMMAND ${CMAKE_COMMAND} ${STORM_3RDPARTY_SOURCE_DIR}/carl "-DSTORM_3RDPARTY_BINARY_DIR=${STORM_3RDPARTY_BINARY_DIR}" "-DBoost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}" "-DBoost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}" "-DUse_CLN_NUMBERS=${SHIPPED_CARL_USE_CLN_NUMBERS}" "-DUse_GINAC=${SHIPPED_CARL_USE_GINAC}" "-DCARL_PORTABLE=${STORM_PORTABLE}"
	WORKING_DIRECTORY ${STORM_3RDPARTY_BINARY_DIR}/carl_download
	OUTPUT_VARIABLE carlconfig_out
	RESULT_VARIABLE carlconfig_result)

	if(NOT carlconfig_result)
		message(STATUS "${carlconfig_out}")
	endif()
	execute_process(
		COMMAND ${CMAKE_COMMAND} --build . --target carl-config
		WORKING_DIRECTORY ${STORM_3RDPARTY_BINARY_DIR}/carl_download
		OUTPUT_VARIABLE carlconfig_out
		RESULT_VARIABLE carlconfig_result
		)
		if(NOT carlconfig_result)
		message(STATUS "${carlconfig_out}")
	endif()
    message(STATUS "Carl - End of config process")

    message(STATUS "Storm - Using shipped version of carl.")
    ExternalProject_Add(
            carl
            SOURCE_DIR ${STORM_3RDPARTY_BINARY_DIR}/carl
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make lib_carl -j${STORM_RESOURCES_BUILD_JOBCOUNT}
            INSTALL_COMMAND make install -j${STORM_RESOURCES_BUILD_JOBCOUNT}
            LOG_BUILD ON
			LOG_INSTALL ON
            BUILD_BYPRODUCTS ${STORM_3RDPARTY_BINARY_DIR}/carl/lib/libcarl${DYNAMIC_EXT}
    )
    include(${STORM_3RDPARTY_BINARY_DIR}/carl/carlConfig.cmake)

    set(STORM_HAVE_CLN ${CARL_USE_CLN_NUMBERS})
    set(STORM_HAVE_GINAC ${CARL_USE_GINAC})

	add_dependencies(resources carl)
    set(carl_INCLUDE_DIR "${STORM_3RDPARTY_BINARY_DIR}/carl/include/")
    set(carl_DIR "${STORM_3RDPARTY_BINARY_DIR}/carl/")
    set(carl_LIBRARIES ${STORM_3RDPARTY_BINARY_DIR}/carl/lib/libcarl${DYNAMIC_EXT})
    set(STORM_HAVE_CARL ON)

    message(STATUS "Storm - Linking with shipped carl ${carl_VERSION} (include: ${carl_INCLUDE_DIR}, library ${carl_LIBRARIES}, CARL_USE_CLN_NUMBERS: ${CARL_USE_CLN_NUMBERS}, CARL_USE_GINAC: ${CARL_USE_GINAC}).")

    # install the carl dynamic library if we built it
    if(MACOSX)
        install(FILES ${STORM_3RDPARTY_BINARY_DIR}/carl/lib/libcarl.${carl_VERSION}${DYNAMIC_EXT} DESTINATION lib)
    else()
        install(FILES ${STORM_3RDPARTY_BINARY_DIR}/carl/lib/libcarl${DYNAMIC_EXT}.${carl_VERSION} DESTINATION lib)
    endif()
endif()

if("${carl_VERSION_MAJOR}.${carl_VERSION_MINOR}" VERSION_EQUAL "14.22")
    # This version is too old for forward declarations and updating requires moving the git,
    # so we warn users but start warning them now.
    set(STORM_CARL_SUPPORTS_FWD_DECL OFF)
    message(WARNING "Uses an outdated repo for Carl. Carl is now hosted at https://github.com/moves-rwth/carl-storm")
elseif("${carl_VERSION_MAJOR}.${carl_VERSION_MINOR}" VERSION_EQUAL "14.23")
    # This version is too old for forward declarations, but we keep supporting it for the moment.
    set(STORM_CARL_SUPPORTS_FWD_DECL OFF)
else()
    set(STORM_CARL_SUPPORTS_FWD_DECL ON)
endif()

if(STORM_USE_CLN_RF AND NOT STORM_HAVE_CLN)
	message(FATAL_ERROR "Cannot use CLN numbers if carl is build without.")
endif()
if(STORM_USE_CLN_RF AND NOT STORM_HAVE_GINAC)
    message(FATAL_ERROR "Cannot use CLN numbers if carl is build without ginac.")
endif()

# The library that needs symbols must be first, then the library that resolves the symbol.
list(APPEND STORM_DEP_IMP_TARGETS lib_carl)
if(STORM_USE_CLN_EA OR STORM_USE_CLN_RF)
    list(APPEND STORM_DEP_IMP_TARGETS GINAC_SHARED CLN_SHARED)
endif()
list(APPEND STORM_DEP_IMP_TARGETS GMPXX_SHARED GMP_SHARED)


#############################################################
##
##      SMT-RAT
##
#############################################################

set(STORM_HAVE_SMTRAT OFF)
if(USE_SMTRAT)
	find_package(smtrat QUIET REQUIRED)
	if(smtrat_FOUND)
		set(STORM_HAVE_SMTRAT ON)
		message(STATUS "Storm - Linking with smtrat.")
    	include_directories("${smtrat_INCLUDE_DIR}")
    	list(APPEND STORM_LINK_LIBRARIES ${smtrat_LIBRARIES})
	else()
		message(FATAL_ERROR "Storm - SMT-RAT was requested but not found")
	endif()
endif()


#############################################################
##
##	gmp
##
#############################################################

get_target_property(GMPXX_LIB GMPXX_SHARED IMPORTED_LIB_LOCATION)
get_target_property(GMP_LIB GMP_SHARED IMPORTED_LIB_LOCATION)
get_target_property(GMPXX_INCLUDE_DIR GMPXX_SHARED INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(GMP_INCLUDE_DIR GMP_SHARED INTERFACE_INCLUDE_DIRECTORIES)
get_filename_component(GMP_LIB_LOCATION ${GMP_LIB} DIRECTORY)
get_filename_component(GMPXX_LIB_LOCATION ${GMPXX_LIB} DIRECTORY)


#############################################################
##
##	cln
##
#############################################################

if(STORM_HAVE_CLN)
	get_target_property(CLN_INCLUDE_DIR CLN_SHARED INTERFACE_INCLUDE_DIRECTORIES)
endif()

#############################################################
##
##	MathSAT (optional)
##
#############################################################

if ("${MSAT_ROOT}" STREQUAL "")
    set(ENABLE_MSAT OFF)
else()
    set(ENABLE_MSAT ON)
endif()

# MathSAT Defines
set(STORM_HAVE_MSAT ${ENABLE_MSAT})
if (ENABLE_MSAT)
    message (STATUS "Storm - Linking with MathSAT.")
    if(${OPERATING_SYSTEM} MATCHES "Linux")
        find_library(MSAT_LIB mathsat PATHS "${MSAT_ROOT}/lib")
        add_imported_library(msat SHARED ${MSAT_LIB} "${MSAT_ROOT}/include")
        list(APPEND STORM_DEP_TARGETS msat_SHARED)
    else()
		# on macOS, the .dylib file has some hard coded path (Version 5.5.4).
		# we thus link statically
        find_library(MSAT_LIB NAMES libmathsat${STATIC_EXT} mathsat PATHS "${MSAT_ROOT}/lib")
        add_imported_library(msat STATIC ${MSAT_LIB} "${MSAT_ROOT}/include")
        list(APPEND STORM_DEP_TARGETS msat_STATIC)
    endif()
endif(ENABLE_MSAT)

#############################################################
##
##	QVBS (Quantitative verification benchmark set)
##
#############################################################

set(STORM_HAVE_QVBS OFF)
if (STORM_LOAD_QVBS)
    ExternalProject_Add(
            download_qvbs
            GIT_REPOSITORY https://github.com/ahartmanns/qcomp.git
            GIT_SHALLOW 1
            PREFIX ${PROJECT_BINARY_DIR}/qvbs
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            LOG_INSTALL ON
    )
    add_dependencies(resources download_qvbs)
    set(STORM_HAVE_QVBS ON)
    ExternalProject_Get_Property(download_qvbs source_dir)
    set(STORM_QVBS_ROOT "${source_dir}/benchmarks")
elseif (NOT "${STORM_QVBS_ROOT}" STREQUAL "")
    if( EXISTS "${STORM_QVBS_ROOT}/index.json")
        set(STORM_HAVE_QVBS ON)
    else()
        message(FATAL_ERROR "No file 'index.json' in provided QVBS root directory ${STORM_QVBS_ROOT}")
    endif()
endif()
if (STORM_HAVE_QVBS)
    message (STATUS "Storm - Enabled inputs from QVBS located at ${STORM_QVBS_ROOT}")
endif()


#############################################################
##
##	Xerces
##
#############################################################

include(${STORM_3RDPARTY_SOURCE_DIR}/include_xerces.cmake)

#############################################################
##
##	SoPlex
##
#############################################################

include(${STORM_3RDPARTY_SOURCE_DIR}/include_soplex.cmake)

#############################################################
##
##	Spot
##
#############################################################

include(${STORM_3RDPARTY_SOURCE_DIR}/include_spot.cmake)

#############################################################
##
##	Sylvan
##
#############################################################

if(STORM_SHIPPED_CARL)
	set(sylvan_dep carl)
else()
	set(sylvan_dep lib_carl)
endif()

if (STORM_DEBUG_SYLVAN)
	set(SYLVAN_BUILD_TYPE "Debug")
	message(WARNING "Storm - Building Sylvan in DEBUG mode.")
else()
	set(SYLVAN_BUILD_TYPE "Release")
endif()

ExternalProject_Add(
        sylvan
        DOWNLOAD_COMMAND ""
        PREFIX "sylvan"
        SOURCE_DIR ${STORM_3RDPARTY_SOURCE_DIR}/sylvan
        CMAKE_ARGS -DPROJECT_NAME=storm -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DSYLVAN_BUILD_DOCS=OFF -DSYLVAN_BUILD_EXAMPLES=OFF -DSYLVAN_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=${SYLVAN_BUILD_TYPE} -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DSYLVAN_GMP=ON -DUSE_CARL=ON -Dcarl_DIR=${carl_DIR} -DBUILD_SHARED_LIBS=OFF
        BINARY_DIR ${STORM_3RDPARTY_BINARY_DIR}/sylvan
        BUILD_IN_SOURCE 0
        INSTALL_COMMAND ""
        INSTALL_DIR ${STORM_3RDPARTY_BINARY_DIR}/sylvan
        LOG_CONFIGURE ON
        LOG_BUILD ON
        DEPENDS ${sylvan_dep}
        BUILD_BYPRODUCTS ${STORM_3RDPARTY_BINARY_DIR}/sylvan/src/libsylvan${STATIC_EXT}
        BUILD_ALWAYS 1
)
# BUILD ALWAYS ON due to: https://stackoverflow.com/questions/46708124/cmake-doesnt-rebuild-externalproject-on-changes

ExternalProject_Get_Property(sylvan source_dir)
ExternalProject_Get_Property(sylvan binary_dir)
set(sylvan_INCLUDE_DIR "${source_dir}/src")
set(sylvan_LIBRARY "${binary_dir}/src/libsylvan${STATIC_EXT}")
message(STATUS "Storm - Using shipped version of sylvan.")
message(STATUS "Storm - Linking with sylvan.")
add_imported_library(sylvan STATIC ${sylvan_LIBRARY} ${sylvan_INCLUDE_DIR})
add_dependencies(sylvan_STATIC sylvan)

list(APPEND STORM_DEP_TARGETS sylvan_STATIC)

#############################################################
##
##	Google Test gtest
##
#############################################################
set(GOOGLETEST_LIB_DIR ${STORM_3RDPARTY_BINARY_DIR}/googletest)
ExternalProject_Add(
        googletest
        #For downloads (may be useful later!)
        #SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
        #TIMEOUT 10
        DOWNLOAD_COMMAND ""
        SOURCE_DIR "${STORM_3RDPARTY_SOURCE_DIR}/googletest"
        # Force the same output paths for debug and release builds so that
        # we know in which place the binaries end up when using the Xcode generator
        CMAKE_ARGS  -Dgtest_force_shared_crt=ON -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${GOOGLETEST_LIB_DIR} -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${GOOGLETEST_LIB_DIR}
        # Disable install step
        INSTALL_COMMAND ""
        BINARY_DIR ${GOOGLETEST_LIB_DIR}
        INSTALL_DIR ${GOOGLETEST_LIB_DIR}
        # Wrap download, configure and build steps in a script to log output
        LOG_CONFIGURE ON
        LOG_BUILD ON
        BUILD_BYPRODUCTS ${GOOGLETEST_LIB_DIR}/lib/libgtest${STATIC_EXT} ${GOOGLETEST_LIB_DIR}/lib/libgtest_main${STATIC_EXT}
)

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIR ${source_dir}/googletest/include)
set(GMOCK_INCLUDE_DIR ${source_dir}/googlemock/include)
# Specify MainTest's link libraries
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBRARIES ${binary_dir}/lib/libgtest${STATIC_EXT} ${binary_dir}/lib/libgtest_main${STATIC_EXT})

add_dependencies(test-resources googletest)
list(APPEND STORM_TEST_LINK_LIBRARIES ${GTEST_LIBRARIES})

#############################################################
##
##	Intel Threading Building Blocks (optional)
##
#############################################################

set(STORM_HAVE_INTELTBB OFF)
if (STORM_USE_INTELTBB)
    find_package(TBB QUIET REQUIRED)

    if (TBB_FOUND)
        message(STATUS "Storm - Found Intel TBB with interface version ${TBB_INTERFACE_VERSION}.")
        message(STATUS "Storm - Linking with Intel TBB in ${TBB_LIBRARY_DIRS}.")
        set(STORM_HAVE_INTELTBB ON)
	# Under Linux libtbb.so contains a linker script to libtbb.so.2 instead of a symlink.
	# This breaks CMake.
        # As a workaround we manually try to add the necessary suffix ".2" to the so file and hope for the best.
	if (LINUX)
		# Check if the library is a linker script which only links to the correct library
		# Read first bytes of file
		file(READ "${TBB_LIBRARY}" TMPTXT LIMIT 128)
		# Check if file starts with "INPUT (libtbb.so.2)"
		if("${TMPTXT}" MATCHES "INPUT \\(libtbb\\.so\\.(.*)\\)")
			# Manually set library by adding the suffix from the linker script.
			# CMAKE_MATCH_1 contains the parsed suffix.
			set(TBB_LIB_LINUX "${TBB_LIBRARY}.${CMAKE_MATCH_1}")
			set(TBB_MALLOC_LIB_LINUX "${TBB_MALLOC_LIBRARY}.${CMAKE_MATCH_1}")
			if (EXISTS "${TBB_LIB_LINUX}")
				set(TBB_LIBRARY ${TBB_LIB_LINUX})
				message(STATUS "Storm - Using Intel TBB library in manually set path ${TBB_LIBRARY}.")
			endif()
			if (EXISTS "${TBB_MALLOC_LIB_LINUX}")
				set(TBB_MALLOC_LIBRARY ${TBB_MALLOC_LIB_LINUX})
				message(STATUS "Storm - Using Intel TBB malloc library in manually set path ${TBB_MALLOC_LIBRARY}.")
			endif()
		endif()
	endif()

        add_imported_library(tbb SHARED ${TBB_LIBRARY} ${TBB_INCLUDE_DIRS})
        list(APPEND STORM_DEP_TARGETS tbb_SHARED)
        add_imported_library(tbb_malloc SHARED ${TBB_MALLOC_LIBRARY} ${TBB_INCLUDE_DIRS})
        list(APPEND STORM_DEP_TARGETS tbb_malloc_SHARED)

    else(TBB_FOUND)
        message(FATAL_ERROR "Storm - TBB was requested, but not found.")
    endif(TBB_FOUND)
endif(STORM_USE_INTELTBB)

#############################################################
##
##	Threads
##
#############################################################

find_package(Threads QUIET REQUIRED)
if (NOT Threads_FOUND)
    message(FATAL_ERROR "Storm - Threads was requested, but not found.")
endif()
include_directories(${THREADS_INCLUDE_DIRS})
list(APPEND STORM_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
if (STORM_USE_COTIRE)
    target_link_libraries(storm_unity ${CMAKE_THREAD_LIBS_INIT})
endif(STORM_USE_COTIRE)


add_custom_target(copy_resources_headers DEPENDS ${PHMAP_BINDIR_HEADERS})
