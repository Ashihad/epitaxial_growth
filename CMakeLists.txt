cmake_minimum_required(VERSION 3.22)

project(
    "EpitaxialGrowth"
    VERSION 1.0.0
    LANGUAGES C CXX)

# Global CMake variables are set here
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(ENABLE_WARNINGS "Enable to add warnings to a target." ON)
option(ENABLE_WARNINGS_AS_ERRORS "Enable to treat warnings as errors." OFF)

option(ENABLE_CLANG_TIDY "Enable to add clang tidy." ON)

option(ENABLE_SANITIZE_ADDR "Enable address sanitize." OFF)
option(ENABLE_SANITIZE_UNDEF "Enable undefined sanitize." OFF)
option(ENABLE_SANITIZE_LEAK "Enable leak sanitize." OFF)
option(ENABLE_SANITIZE_THREAD "Enable thread sanitize." OFF)

option(ENABLE_CLANG_FORMAT "Enable to add clang-format." ON)
option(ENABLE_CLANG_TIDY "Enable to add clang-tidy." ON)

# Project/Library Names
set(EXECUTABLE_NAME "main")

# CMAKE MODULES
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
include(ConfigSafeGuards)
include(AddGitSubmodule)
include(Docs)
include(Tools)

if(ENABLE_WARNINGS)
    include(Warnings)
endif()

add_clang_format_target()

if(ENABLE_SANITIZE_ADDR
   OR ENABLE_SANITIZE_UNDEF
   OR ENABLE_SANITIZE_LEAK
   OR ENABLE_SANITIZE_THREAD)
    include(Sanitizer)
    add_sanitizer_flags()
endif()

# EXTERNAL LIBRARIES
message(STATUS "Using CPM")
include(CPM)
cpmaddpackage("gl:libeigen/eigen#3.4")

# enable LTO manually
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

# BUILDING
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

set(MAIN_FILE ${CMAKE_SOURCE_DIR}/src/main.cpp)
include_directories(${CMAKE_SOURCE_DIR}/hdr)

file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM SRC_FILES ${MAIN_FILE})

message(STATUS "Eigen include directory: ${eigen_INCLUDE_DIRS}")

# Build all nonmain targets
add_library(objects OBJECT ${SRC_FILES})
if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        objects
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()
# Disable warnings specifically for Eigen by including it as system library
target_include_directories(objects SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/build/_deps/eigen-src/)
if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target(objects)
endif()

# Build main
add_executable(epitaxal_mc $<TARGET_OBJECTS:objects> ${MAIN_FILE})
if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        epitaxal_mc
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()
# Disable warnings specifically for Eigen by including it as system library
target_include_directories(epitaxal_mc SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/build/_deps/eigen-src/)
if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target(epitaxal_mc)
endif()
