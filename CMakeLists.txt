cmake_minimum_required(VERSION 3.22)

project(
  "EpitaxialGrowth"
  VERSION 1.0.0
  LANGUAGES C CXX
)

# Global CMake variables are set here
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Options
option(ENABLE_WARNINGS "Enable to add warnings to a target." ON)
option(ENABLE_WARNINGS_AS_ERRORS "Enable to treat warnings as errors." OFF)

option(ENABLE_TESTING "Enable a Unit Testing build." OFF)
option(ENABLE_COVERAGE "Enable a Code Coverage build." OFF)

option(ENABLE_SANITIZE_ADDR "Enable address sanitize." OFF)
option(ENABLE_SANITIZE_UNDEF "Enable undefined sanitize." OFF)
option(ENABLE_SANITIZE_LEAK "Enable leak sanitize." OFF)
option(ENABLE_SANITIZE_THREAD "Enable thread sanitize." OFF)

option(ENABLE_LTO "Enable to add Link Time Optimization." ON)

# CMAKE MODULES
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
include(ConfigSafeGuards)
include(AddGitSubmodule)
include(Docs)

if(ENABLE_WARNINGS)
    include(Warnings)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(ENABLE_SANITIZE_ADDR ON)
  set(ENABLE_SANITIZE_UNDEF ON)
  set(ENABLE_SANITIZE_LEAK OFF)
  set(ENABLE_SANITIZE_THREAD OFF)
  include(Sanitizer)
  add_sanitizer_flags()

  set(ENABLE_TESTING ON)
  set(ENABLE_COVERAGE ON)
  include(CodeCoverage)
  append_coverage_compiler_flags()

  include(CTest)
  add_subdirectory(tests)
else()
  set(ENABLE_SANITIZE_ADDR OFF)
  set(ENABLE_SANITIZE_UNDEF OFF)
  set(ENABLE_SANITIZE_LEAK OFF)
  set(ENABLE_SANITIZE_THREAD OFF)

  set(ENABLE_TESTING OFF)
  set(ENABLE_COVERAGE OFF)

  include(PerformanceFlags)
  append_performance_compiler_flags()
endif()

# Project/Library Names
set(EXECUTABLE_NAME "main")

if(ENABLE_LTO)
  include(LTO)
endif()

# EXTERNAL LIBRARIES
message(STATUS "Using CPM")
include(CPM)
# cpmaddpackage("gl:libeigen/eigen#3.4")
if(ENABLE_TESTING)
  cpmaddpackage("gh:google/googletest#v1.15.x")
endif()

# BUILDING
add_subdirectory(src)
add_subdirectory(app)

install(
  TARGETS ${EXECUTABLE_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(
  TARGETS "atom" "atom_diffused" "file" "linalg" "rng" "sim"
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
