set(APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
set(
  APP_INCLUDES
  "${CMAKE_SOURCE_DIR}/src/atom"
  "${CMAKE_SOURCE_DIR}/src/atom_diffused"
  "${CMAKE_SOURCE_DIR}/src/config"
  "${CMAKE_SOURCE_DIR}/src/file"
  "${CMAKE_SOURCE_DIR}/src/linalg"
  "${CMAKE_SOURCE_DIR}/src/physics"
  "${CMAKE_SOURCE_DIR}/src/rng"
  "${CMAKE_SOURCE_DIR}/src/sim"
  "${CMAKE_SOURCE_DIR}/src/utils"
)


add_executable(${EXECUTABLE_NAME} ${APP_SOURCES})
target_include_directories(${EXECUTABLE_NAME} PUBLIC ${APP_INCLUDES})

target_link_libraries(
  ${EXECUTABLE_NAME}
  PRIVATE
  "atom"
  "atom_diffused"
  "file"
  "linalg"
  "rng"
  "sim"
)

if(${ENABLE_WARNINGS})
  target_set_warnings(
    TARGET
    ${EXECUTABLE_NAME}
    ENABLE
    ${ENABLE_WARNINGS}
    AS_ERRORS
    ${ENABLE_WARNINGS_AS_ERRORS}
  )
endif()

if(${ENABLE_LTO})
  target_enable_lto(
    TARGET
    ${EXECUTABLE_NAME}
    ENABLE
    ON
  )
endif()
