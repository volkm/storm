set(STORM_HAVE_GMM OFF)
set(STORM_SHIPPED_GMM OFF)
set(STORM_GMM_NEEDS_BLAS OFF)
set(STORM_GMM_NEEDS_LAPACK OFF)

if(NOT STORM_DISABLE_GMM)
    # Try to find gmm on the system
    if (NOT STORM_GMM_FORCE_SHIPPED)
        if (NOT "${GMM_ROOT}" STREQUAL "")
            message(STATUS "Storm - Searching for gmm in ${GMM_ROOT}")
            find_path(GMM_INCLUDE_DIR gmm/gmm_kernel.h PATHS ${GMM_ROOT} NO_DEFAULT_PATH)
            if(GMM_INCLUDE_DIR MATCHES "-NOTFOUND")
                unset(GMM_INCLUDE_DIR)
            endif()
        endif()
        if (NOT GMM_INCLUDE_DIR)
            find_path(GMM_INCLUDE_DIR gmm/gmm_kernel.h)
        endif()

        if (GMM_INCLUDE_DIR)
            set(_GMM_VERSION_STRING "unknown version")
            file(STRINGS "${GMM_INCLUDE_DIR}/gmm/gmm_arch_config.h" _GMM_VERSION_LINE REGEX "[ \t]*#define GMM_VERSION \"")
            if(_GMM_VERSION_LINE)
                string(REGEX REPLACE "^[ \t]*#define GMM_VERSION \"([^\"]*)\"" "\\1" _GMM_VERSION_STRING "${_GMM_VERSION_LINE}")
            endif()

            # The shipped gmm is built without BLAS/LAPACK support (see gmm55.patch).
            # A system gmm might be built with GMM_USES_BLAS/GMM_USES_LAPACK defined,
            # in which case BLAS/LAPACK have to be linked against.
            file(STRINGS "${GMM_INCLUDE_DIR}/gmm/gmm_arch_config.h" _GMM_USES_BLAS REGEX "^#define GMM_USES_BLAS")
            file(STRINGS "${GMM_INCLUDE_DIR}/gmm/gmm_arch_config.h" _GMM_USES_LAPACK REGEX "^#define GMM_USES_LAPACK")
            if(_GMM_USES_BLAS)
                find_package(BLAS QUIET)
                if(NOT BLAS_FOUND)
                    message(FATAL_ERROR "Storm - The system gmm requires linking with BLAS, but no BLAS library was found. Install a BLAS library or force the shipped gmm version with '-DSTORM_GMM_FORCE_SHIPPED=ON'.")
                endif()
                set(STORM_GMM_NEEDS_BLAS ON)
            endif()
            if(_GMM_USES_LAPACK)
                find_package(LAPACK QUIET)
                if(NOT LAPACK_FOUND)
                    message(FATAL_ERROR "Storm - The system gmm requires linking with LAPACK, but no LAPACK library was found. Install a LAPACK library or force the shipped gmm version with '-DSTORM_GMM_FORCE_SHIPPED=ON'.")
                endif()
                set(STORM_GMM_NEEDS_LAPACK ON)
            endif()

            add_library(gmm INTERFACE IMPORTED)
            if(STORM_GMM_NEEDS_BLAS OR STORM_GMM_NEEDS_LAPACK)
                set_target_properties(gmm PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${GMM_INCLUDE_DIR}"
                        INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}"
                )
            else()
                set_target_properties(gmm PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GMM_INCLUDE_DIR}")
            endif()

            message(STATUS "Storm - Using system version of gmm ${_GMM_VERSION_STRING} (include: ${GMM_INCLUDE_DIR}, linking BLAS: ${STORM_GMM_NEEDS_BLAS}, LAPACK: ${STORM_GMM_NEEDS_LAPACK}).")
            list(APPEND STORM_DEP_IMP_TARGETS gmm)

            set(STORM_HAVE_GMM ON)
        endif()
    endif()

    if(NOT STORM_HAVE_GMM)
        set(GMM_SHIPPED_VERSION "5.5")
        message(STATUS "Storm - Including shipped gmm ${GMM_SHIPPED_VERSION}.")
        ExternalProject_Add(
                gmm_src
                URL https://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-${GMM_SHIPPED_VERSION}.tar.gz
                SOURCE_DIR ${STORM_3RDPARTY_BINARY_DIR}/gmm
                # Try to apply patch and otherwise check that patch was already applied
                # Needed because while -N ignores already applied patches it still returns with exit code 1
                PATCH_COMMAND patch -N -p1 -s -i ${STORM_3RDPARTY_SOURCE_DIR}/patches/gmm55.patch || patch --dry-run -N -p1 -s -i ${STORM_3RDPARTY_SOURCE_DIR}/patches/gmm55.patch | grep "Reversed (or previously applied) patch detected!"
                UPDATE_COMMAND ""
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND ""
                LOG_DOWNLOAD ON
                LOG_INSTALL ON
                LOG_OUTPUT_ON_FAILURE ON
        )
        add_library(gmm INTERFACE) # Not imported, we are in control of the sources.
        add_dependencies(gmm gmm_src)
        add_dependencies(storm_resources gmm)
        target_include_directories(gmm INTERFACE
                $<BUILD_INTERFACE:${STORM_3RDPARTY_BINARY_DIR}/gmm/include>
                $<INSTALL_INTERFACE:${STORM_RESOURCE_INCLUDE_INSTALL_DIR}>
        )
        install(TARGETS gmm EXPORT storm_Targets)
        install(DIRECTORY ${STORM_3RDPARTY_BINARY_DIR}/gmm/include/gmm
                DESTINATION ${STORM_RESOURCE_INCLUDE_INSTALL_DIR}
                FILES_MATCHING PATTERN "*.h" PATTERN ".git" EXCLUDE)
        list(APPEND STORM_DEP_TARGETS gmm)
        set(STORM_HAVE_GMM ON)
        set(STORM_SHIPPED_GMM ON)

        message(STATUS "Storm - Using shipped version of gmm ${GMM_SHIPPED_VERSION}.")
    endif()
else()
    message(STATUS "Storm - Not linking with gmm.")
    set(STORM_HAVE_GMM OFF)
endif()
