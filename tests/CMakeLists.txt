include(GoogleTest)
enable_testing()

add_executable("main_UT" "test_linalg.cpp" "test_physics.cpp" "test_rng.cpp" "test_sim.cpp")
target_link_libraries("main_UT" PUBLIC "linalg" "physics" "rng" "sim")
target_link_libraries("main_UT" PRIVATE GTest::gtest_main)
if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "main_UT"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()
gtest_discover_tests("main_UT")

if(ENABLE_COVERAGE)
    set(COVERAGE_MAIN "coverage")
    set(COVERAGE_EXCLUDES
        "${PROJECT_SOURCE_DIR}/app/*"
        "${PROJECT_SOURCE_DIR}/cmake/*"
        "${PROJECT_SOURCE_DIR}/docs/*"
        "${PROJECT_SOURCE_DIR}/external/*"
        "${PROJECT_SOURCE_DIR}/tests/*"
        "${PROJECT_SOURCE_DIR}/build/*"
        "/usr/include/*")
    set(COVERAGE_EXTRA_FLAGS)
    set(COVERAGE_DEPENDENCIES "main_UT")

    setup_target_for_coverage_gcovr_html(
        NAME
        ${COVERAGE_MAIN}
        EXECUTABLE
        ctest
        DEPENDENCIES
        ${COVERAGE_DEPENDENCIES}
        BASE_DIRECTORY
        ${CMAKE_SOURCE_DIR})
endif()