include(CodeCoverage)
append_coverage_compiler_flags()


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

add_executable(test_chyps_serial)
target_include_directories(test_chyps_serial PRIVATE "${PROJECT_SOURCE_DIR}")
target_link_libraries(test_chyps_serial chyps gtest gmock)

set_target_properties(test_chyps_serial
    PROPERTIES
    CXX_STANDARD 17
)

add_executable(test_chyps_parallel)
target_include_directories(test_chyps_parallel PRIVATE "${PROJECT_SOURCE_DIR}")
target_include_directories(test_chyps_parallel PRIVATE "${GTEST_MPI_SINGLE_HEADER}")
target_link_libraries(test_chyps_parallel chyps gtest gmock)

set_target_properties(test_chyps_parallel
    PROPERTIES
    CXX_STANDARD 17
)

add_executable(test_chyps_verification)
target_include_directories(test_chyps_verification PRIVATE "${PROJECT_SOURCE_DIR}")
target_include_directories(test_chyps_verification PRIVATE "${GTEST_MPI_SINGLE_HEADER}")
target_link_libraries(test_chyps_verification chyps gtest gmock)

set_target_properties(test_chyps_verification
    PROPERTIES
    CXX_STANDARD 17
)



include(GoogleTest)

add_subdirectory("./unit")
add_subdirectory("./verification")

setup_target_for_coverage_gcovr_html(
NAME test_chyps_coverage
EXECUTABLE tests/test_chyps_serial
DEPENDENCIES test_chyps_serial test_chyps_parallel
)

add_custom_command(TARGET test_chyps_serial POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${CMAKE_CURRENT_SOURCE_DIR}/unit/serial/data" "${CMAKE_CURRENT_BINARY_DIR}/unit/serial/data")
add_custom_command(TARGET test_chyps_parallel POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${CMAKE_CURRENT_SOURCE_DIR}/unit/parallel/data" "${CMAKE_CURRENT_BINARY_DIR}/unit/parallel/data")
add_custom_command(TARGET test_chyps_verification POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${CMAKE_CURRENT_SOURCE_DIR}/verification/data" "${CMAKE_CURRENT_BINARY_DIR}/verification/data")

gtest_add_tests(TARGET test_chyps_serial WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
gtest_add_tests(TARGET test_chyps_parallel WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" TEST_SUFFIX "_1core" TEST_LIST core1)
gtest_add_tests(TARGET test_chyps_parallel WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" TEST_SUFFIX "_2core" TEST_LIST core2)
gtest_add_tests(TARGET test_chyps_parallel WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" TEST_SUFFIX "_4core" TEST_LIST core4)
set_tests_properties(${core1} PROPERTIES PROCESSORS 1)
set_tests_properties(${core2} PROPERTIES PROCESSORS 1)
set_tests_properties(${core4} PROPERTIES PROCESSORS 1)

