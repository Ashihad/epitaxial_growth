include(GoogleTest)
enable_testing()

set(UT_EXEC "main_UT")

add_executable(
  ${UT_EXEC}
  "test_linalg.cpp"
  "test_physics.cpp"
  "test_rng.cpp"
  "test_sim.cpp"
)

target_include_directories(${UT_EXEC} PUBLIC .)

target_link_libraries(
  ${UT_EXEC}
  PUBLIC
  "atom"
  "atom_diffused"
  "file"
  "linalg"
  "rng"
  "sim"
)
target_link_libraries(${UT_EXEC} PRIVATE GTest::gtest_main)

if(${ENABLE_WARNINGS})
  target_set_warnings(
    TARGET
    ${UT_EXEC}
    ENABLE
    ${ENABLE_WARNINGS}
    AS_ERRORS
    ${ENABLE_WARNINGS_AS_ERRORS}
  )
endif()
gtest_discover_tests(${UT_EXEC})

if(ENABLE_COVERAGE)
  set(COVERAGE_MAIN "coverage")
  set(COVERAGE_EXCLUDES
    "${PROJECT_SOURCE_DIR}/app/*"
    "${PROJECT_SOURCE_DIR}/cmake/*"
    "${PROJECT_SOURCE_DIR}/docs/*"
    "${PROJECT_SOURCE_DIR}/external/*"
    "${PROJECT_SOURCE_DIR}/tests/*"
    "${PROJECT_SOURCE_DIR}/build/*"
    "${PROJECT_SOURCE_DIR}/debug/*"
    "${PROJECT_SOURCE_DIR}/release/*"
    "/usr/include/*"
  )
  set(COVERAGE_EXTRA_FLAGS)
  set(COVERAGE_DEPENDENCIES ${UT_EXEC})

  setup_target_for_coverage_gcovr_html(
    NAME
    ${COVERAGE_MAIN}
    EXECUTABLE
    ctest
    DEPENDENCIES
    ${COVERAGE_DEPENDENCIES}
    BASE_DIRECTORY
    ${CMAKE_SOURCE_DIR}
  )
endif()