set(LIBRARY_NAME "file")
set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/FileHandler.cpp")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/FileHandler.hpp")
set(
  LIBRARY_INCLUDES
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

add_library(${LIBRARY_NAME} STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories(${LIBRARY_NAME} PUBLIC ${LIBRARY_INCLUDES})
# target_link_libraries(${LIBRARY_NAME} PRIVATE "physics" "rng")

if(${ENABLE_WARNINGS})
  target_set_warnings(
    TARGET
    ${LIBRARY_NAME}
    ENABLE
    ${ENABLE_WARNINGS}
    AS_ERRORS
    ${ENABLE_WARNINGS_AS_ERRORS}
  )
endif()

if(${ENABLE_LTO})
  target_enable_lto(
    TARGET
    ${LIBRARY_NAME}
    ENABLE
    ON
  )
endif()
