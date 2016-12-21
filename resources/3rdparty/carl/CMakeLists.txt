project(carlext)
cmake_minimum_required(VERSION 3.3)
include(ExternalProject)

option(STORM_3RDPARTY_BINARY_DIR "3rd party bin dir")

message(STORM_3RDPARTY_BINARY_DIR: ${STORM_3RDPARTY_BINARY_DIR})

ExternalProject_Add(carl-config
	GIT_REPOSITORY https://github.com/smtrat/carl
        GIT_TAG master
	PREFIX here
	SOURCE_DIR source_dir
	BINARY_DIR ${STORM_3RDPARTY_BINARY_DIR}/carl
	CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS} -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS} -DBoost_NO_SYSTEM_PATHS=1 -DEXPORT_TO_CMAKE=0 -DUSE_CLN_NUMBERS=1 -DCMAKE_SKIP_INSTALL_ALL_DEPENDENCY=1 -DCMAKE_INSTALL_PREFIX:PATH=${STORM_3RDPARTY_BINARY_DIR}/carl 
        BUILD_IN_SOURCE 0
        LOG_UPDATE OFF
        LOG_CONFIGURE OFF
        BUILD_COMMAND "" # Disable build step.
	INSTALL_COMMAND "" # Disable install step too.
)
add_custom_target(build-carl)
add_dependencies(build-carl carl-config)

message("done")