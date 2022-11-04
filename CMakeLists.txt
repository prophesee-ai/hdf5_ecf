cmake_minimum_required(VERSION 3.14.0)
cmake_policy(SET CMP0087 NEW)
project(hdf5_ecf VERSION 1.0.0)
set(HDF5_ECF_PROJECT_VERSION "${PROJECT_VERSION}" CACHE INTERNAL "")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release
        RelWithDebInfo MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(HDF5_ECF_PLUGIN_BUILD_PATH "${CMAKE_BINARY_DIR}/lib/hdf5/plugin" CACHE INTERNAL "" FORCE)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT CMAKE_INSTALL_PREFIX)
    if(APPLE)
        set(HDF5_ECF_PLUGIN_INSTALL_PATH "/opt/local/lib/hdf5/" CACHE INTERNAL "" FORCE)
    elseif(WIN32)
        set(HDF5_ECF_PLUGIN_INSTALL_PATH "$ENV{ALLUSERSPROFILE}/hdf5/lib/plugin/" CACHE INTERNAL "" FORCE)
    else()
        set(HDF5_ECF_PLUGIN_INSTALL_PATH "/usr/local/hdf5/lib/plugin" CACHE INTERNAL "" FORCE)
    endif()
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT CMAKE_INSTALL_PREFIX)
    set(HDF5_ECF_PLUGIN_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/lib/hdf5/plugin" CACHE INTERNAL "" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT CMAKE_INSTALL_PREFIX)

if(WIN32)
    set(CMAKE_DEBUG_POSTFIX "_d")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED)

# libhdf5_ecf_codec
add_library(hdf5_ecf_codec SHARED ecf_codec.cpp)

include(CTest)
if (BUILD_TESTING)
    find_package(GTest REQUIRED)
    add_executable(gtest_ecf_codec ecf_codec_test.cpp)
    target_link_libraries(gtest_ecf_codec PRIVATE hdf5_ecf_codec GTest::gtest GTest::gtest_main)
    add_test(
        NAME ecf_codec
        COMMAND $<TARGET_FILE:gtest_ecf_codec> --gtest_color=yes --gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}/gtest/ecf_codec.xml
    )
endif ()

# targets generation and installation
include(GenerateExportHeader)
generate_export_header(hdf5_ecf_codec BASE_NAME ecf_codec)
target_include_directories(hdf5_ecf_codec 
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/..>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
        "$<INSTALL_INTERFACE:include>"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ecf_codec_export.h"
        DESTINATION "include/hdf5_ecf"
        COMPONENT hdf5-ecf-codec-dev
)

# LICENSE
install(FILES "${PROJECT_SOURCE_DIR}/LICENSE"
        DESTINATION "share/hdf5_ecf"
        COMPONENT hdf5-ecf-codec-dev
)

set_target_properties(hdf5_ecf_codec
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
)
install(TARGETS hdf5_ecf_codec
        EXPORT hdf5_ecf-targets
        RUNTIME
            DESTINATION "bin"
            COMPONENT hdf5-ecf-codec-lib
        ARCHIVE
            DESTINATION "lib"
            COMPONENT hdf5-ecf-codec-lib
        LIBRARY
            DESTINATION "lib"
            COMPONENT hdf5-ecf-codec-lib
            NAMELINK_SKIP
)
install(TARGETS hdf5_ecf_codec
        EXPORT hdf5_ecf-targets
        LIBRARY
            DESTINATION "lib"
            COMPONENT hdf5-ecf-codec-dev
            NAMELINK_ONLY
)
install(FILES ecf_codec.h
        DESTINATION "include/hdf5_ecf"
        COMPONENT hdf5-ecf-codec-dev
)
export(EXPORT hdf5_ecf-targets FILE hdf5_ecf-targets.cmake)
install(EXPORT hdf5_ecf-targets
        FILE hdf5_ecf-targets.cmake
        DESTINATION share/cmake/hdf5_ecf
        COMPONENT hdf5-ecf-codec-dev
)

# Create and install configuration files
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "hdf5_ecf-config.cmake.in"
    "hdf5_ecf-config.cmake"
    INSTALL_DESTINATION share/cmake/hdf5_ecf
)
write_basic_package_version_file(
    "hdf5_ecf-config-version.cmake"
    COMPATIBILITY ExactVersion
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hdf5_ecf-config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/hdf5_ecf-config-version.cmake"
        DESTINATION share/cmake/hdf5_ecf
        COMPONENT hdf5-ecf-codec-dev
)

# HDF5 ECF filter 
install(FILES ecf_h5filter.h
        DESTINATION "include/hdf5_ecf"
        COMPONENT hdf5-plugin-ecf-local-install
)

# HDF5 ECF plugin 
add_library(hdf5_ecf_plugin SHARED ecf_h5plugin.cpp)
set_target_properties(hdf5_ecf_plugin PROPERTIES OUTPUT_NAME H5Zecf)
target_link_libraries(hdf5_ecf_plugin PUBLIC hdf5_ecf_codec ${HDF5_LIBRARIES})
target_include_directories(hdf5_ecf_plugin PUBLIC ${HDF5_INCLUDE_DIRS})
file(MAKE_DIRECTORY "${HDF5_ECF_PLUGIN_BUILD_PATH}")
add_custom_command(
    TARGET hdf5_ecf_plugin POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:hdf5_ecf_plugin>" "${HDF5_ECF_PLUGIN_BUILD_PATH}"
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:hdf5_ecf_codec>" "${HDF5_ECF_PLUGIN_BUILD_PATH}"
)

install(TARGETS hdf5_ecf_plugin
        RUNTIME DESTINATION "${HDF5_ECF_PLUGIN_INSTALL_PATH}" COMPONENT hdf5-plugin-ecf-local-install
        LIBRARY DESTINATION "${HDF5_ECF_PLUGIN_INSTALL_PATH}" COMPONENT hdf5-plugin-ecf-local-install
)
install(CODE "message(NOTICE \"\n\nNote:\n=====\nTo use the HDF5 ECF plugin, it may be necessary to append the path '${HDF5_ECF_PLUGIN_INSTALL_PATH}' to the HDF5_PLUGIN_PATH environment variable in your session or permanently.\n\")"
        COMPONENT hdf5-plugin-ecf-local-install
)
set (HDF5_ECF_PLUGIN_DEB_INSTALL_PATH "/usr/lib/${CMAKE_CXX_LIBRARY_ARCHITECTURE}/hdf5/plugins" CACHE INTERNAL "" FORCE)
if (UNIX AND NOT APPLE)
    find_program(_lsb_release_exec lsb_release)
    execute_process(COMMAND ${_lsb_release_exec} -cs
        OUTPUT_VARIABLE _ubuntu_platform
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if ("${_ubuntu_platform}" STREQUAL "jammy")
        set (HDF5_ECF_PLUGIN_DEB_INSTALL_PATH "/usr/lib/${CMAKE_CXX_LIBRARY_ARCHITECTURE}/hdf5/serial/plugins" CACHE INTERNAL "" FORCE)
    endif ()
endif (UNIX AND NOT APPLE)
install(TARGETS hdf5_ecf_plugin
        RUNTIME DESTINATION ${HDF5_ECF_PLUGIN_DEB_INSTALL_PATH} COMPONENT hdf5-plugin-ecf EXCLUDE_FROM_ALL
        LIBRARY DESTINATION ${HDF5_ECF_PLUGIN_DEB_INSTALL_PATH} COMPONENT hdf5-plugin-ecf EXCLUDE_FROM_ALL
)
file(
    WRITE "${CMAKE_CURRENT_BINARY_DIR}/postinst" "\
#!/usr/bin/env sh\n\
echo \"\"\n\
echo \"Note:\"\n\
echo \"=====\"\n\
echo \"To use the HDF5 ECF plugin, it may be necessary to run 'set HDF5_PLUGIN_PATH=${HDF5_ECF_PLUGIN_DEB_INSTALL_PATH}:\\$\{HDF5_PLUGIN_PATH\}' in your session or shell initialization script.\"\n\
echo \"\"\n\
")
file(
    COPY "${CMAKE_CURRENT_BINARY_DIR}/postinst"
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
if (WIN32)
    # on windows, we put the DLL next to the plugin with a symlink to avoid having to mess with PATH
    install(CODE
        "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink \"$<TARGET_FILE:hdf5_ecf_codec>\" \"${HDF5_ECF_PLUGIN_INSTALL_PATH}/$<TARGET_FILE_NAME:hdf5_ecf_codec>\")"
    )
endif()


# CPack
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "")
set(CPACK_PACKAGE_VENDOR "Prophesee")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION_FULL})
set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar") # Variable introduced in cmake 3.7 (see https://cmake.org/cmake/help/v3.8/release/3.7.html#cpack)

# Follow Debian policy on control files' permissions
# cf https://cmake.org/cmake/help/v3.10/module/CPackDeb.html#variable:CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION
set(CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION TRUE)
set(CPACK_DEBIAN_ENABLE_COMPONENT_DEPENDS ON)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "support@prophesee.ai Copyright (c) Prophesee S.A.") #required
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE https://support.prophesee.ai)
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_COMPONENTS_ALL hdf5-ecf-codec-lib hdf5-ecf-codec-dev hdf5-plugin-ecf)

# Library
set(CPACK_COMPONENT_HDF5-ECF-CODEC-LIB_DESCRIPTION "HDF5 ECF Codec library")
set(CPACK_DEBIAN_HDF5-ECF-CODEC-LIB_FILE_NAME "hdf5-ecf-codec-lib-${PROJECT_VERSION}.deb")
set(CPACK_DEBIAN_HDF5-ECF-CODEC-LIB_PACKAGE_NAME "hdf5-ecf-codec-lib")
# Dev
set(CPACK_COMPONENT_HDF5-ECF-CODEC-DEV_DESCRIPTION "HDF5 ECF Codec (C++) files")
set(CPACK_COMPONENT_HDF5-ECF-CODEC-DEV_DEPENDS hdf5-ecf-codec-lib)
set(CPACK_DEBIAN_HDF5-ECF-CODEC-DEV_FILE_NAME "hdf5-ecf-codec-dev-${PROJECT_VERSION}.deb")
set(CPACK_DEBIAN_HDF5-ECF-CODEC-DEV_PACKAGE_NAME "hdf5-ecf-codec-dev")
set(CPACK_DEBIAN_HDF5-ECF-CODEC-DEV_PACKAGE_DEPENDS libhdf5-dev)
# Plugin
set(CPACK_COMPONENT_HDF5-PLUGIN-ECF_DESCRIPTION "HDF5 ECF Plugin")
set(CPACK_DEBIAN_HDF5-PLUGIN-ECF_FILE_NAME "hdf5-plugin-ecf-${PROJECT_VERSION}.deb")
set(CPACK_DEBIAN_HDF5-PLUGIN-ECF_PACKAGE_NAME "hdf5-plugin-ecf")
set(CPACK_COMPONENT_HDF5-PLUGIN-ECF_DEPENDS hdf5-ecf-codec-lib)
set(CPACK_DEBIAN_HDF5-PLUGIN-ECF_PACKAGE_CONTROL_EXTRA "${CPACK_DEBIAN_HDF5-PLUGIN-ECF_PACKAGE_CONTROL_EXTRA};${CMAKE_CURRENT_BINARY_DIR}/postinst;")


include(CPack) # This has to be at the end
