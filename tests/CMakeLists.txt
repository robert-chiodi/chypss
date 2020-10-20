# If path to GoogleTest not given, then download from git
if(NOT GOOGLETEST_DIR)
  # Download and unpack googletest at configure time
  configure_file("${CMAKE_SOURCE_DIR}/cmake/Utils/GoogleTestDownload.cmake" googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                   ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                   EXCLUDE_FROM_ALL)
else()
  add_subdirectory("${GOOGLETEST_DIR}"
                   "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
                   EXCLUDE_FROM_ALL)
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories(SYSTEM "${gtest_SOURCE_DIR}/include")
endif()

add_executable(test_chyps)
target_include_directories(test_chyps PRIVATE "${PROJECT_SOURCE_DIR}")
target_link_libraries(test_chyps chyps gtest gmock gtest_main)
set_target_properties(test_chyps
    PROPERTIES
    CXX_STANDARD 11
)
include(GoogleTest)

add_subdirectory("./unit")

setup_target_for_coverage_gcovr_html(
NAME test_chyps_coverage
EXECUTABLE tests/test_chyps
DEPENDENCIES test_chyps
)
