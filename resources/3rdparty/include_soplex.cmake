if(STORM_USE_SOPLEX)
    find_package(soplex)

    if(SOPLEX_FOUND)
        get_target_property(soplexLOC libsoplex-pic LOCATION)
        get_target_property(soplexINCLUDE libsoplex-pic INTERFACE_INCLUDE_DIRECTORIES)
        MESSAGE(STATUS "Storm - Linking with SoPlex: (library: ${soplexLOC}; include: ${soplexINCLUDE})")
        list(APPEND STORM_DEP_TARGETS libsoplex-pic)
    else()
        message(FATAL_ERROR "SoPlex library requested but was not found.")
    endif()
    set(STORM_HAVE_SOPLEX ${SOPLEX_FOUND})
else()
    set(STORM_HAVE_SOPLEX OFF)
endif()

