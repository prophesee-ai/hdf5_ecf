cmake_minimum_required(VERSION 3.14.0)
cmake_policy(SET CMP0087 NEW)
project(hdf5_ecf VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release
        RelWithDebInfo MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(APPLE)
    set(PLUGIN_INSTALL_PATH "/opt/local/lib/hdf5/" CACHE PATH "HDF5 plugin installation path")
elseif(WIN32)
    set(PLUGIN_INSTALL_PATH "$ENV{ALLUSERSPROFILE}/hdf5/lib/plugin/" CACHE PATH "HDF5 plugin installation path")
else()
    set(PLUGIN_INSTALL_PATH "/usr/local/hdf5/lib/plugin" CACHE PATH "HDF5 plugin installation path")
endif()

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
    target_link_libraries(gtest_ecf_codec PRIVATE hdf5_ecf_codec GTest::GTest GTest::Main)
    add_test(
        NAME ecf_codec
        COMMAND $<TARGET_FILE:gtest_ecf_codec> --gtest_color=yes --gtest_output=xml:${CMAKE_BINARY_DIR}/gtest/ecf_codec.xml
    )
endif ()

# targets generation and installation
include(GenerateExportHeader)
generate_export_header(hdf5_ecf_codec BASE_NAME ecf_codec)
target_include_directories(hdf5_ecf_codec 
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
        "$<INSTALL_INTERFACE:include>"
)
install(FILES "${CMAKE_BINARY_DIR}/ecf_codec_export.h"
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
install(FILES "${CMAKE_BINARY_DIR}/hdf5_ecf-config.cmake" "${CMAKE_BINARY_DIR}/hdf5_ecf-config-version.cmake"
        DESTINATION share/cmake/hdf5_ecf
        COMPONENT hdf5-ecf-codec-dev
)

# HDF5 ECF filter 
install(FILES ecf_h5filter.h
        DESTINATION "include/hdf5_ecf"
        COMPONENT hdf5-plugin-ecf
)

# HDF5 ECF plugin 
add_library(hdf5_ecf_plugin SHARED ecf_h5plugin.cpp)
set_target_properties(hdf5_ecf_plugin PROPERTIES OUTPUT_NAME H5Zecf)
target_link_libraries(hdf5_ecf_plugin PUBLIC hdf5_ecf_codec ${HDF5_LIBRARIES})
target_include_directories(hdf5_ecf_plugin PUBLIC ${HDF5_INCLUDE_DIRS})
install(TARGETS hdf5_ecf_plugin
        RUNTIME DESTINATION ${PLUGIN_INSTALL_PATH} COMPONENT hdf5-plugin-ecf
        LIBRARY DESTINATION ${PLUGIN_INSTALL_PATH} COMPONENT hdf5-plugin-ecf
)
if (WIN32)
    # on windows, we put the DLL next to the plugin with a symlink to avoid having to mess with PATH
    install(CODE
        "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE:hdf5_ecf_codec> ${PLUGIN_INSTALL_PATH}/$<TARGET_FILE_NAME:hdf5_ecf_codec>)"
    )
endif()

