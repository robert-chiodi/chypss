cmake_minimum_required(VERSION 3.9...3.18)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(CHyPS VERSION 0.1
               DESCRIPTION "MFEM-Based Heat-Solver"
               LANGUAGES CXX)


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
"${PROJECT_SOURCE_DIR}/cmake/Utils")

find_package(MPI REQUIRED)
find_package(HYPRE REQUIRED)
find_package(MFEM REQUIRED)
find_package(METIS REQUIRED)
find_package(PRECICE REQUIRED)
find_package(Boost 1.65.1 REQUIRED
  COMPONENTS filesystem log log_setup program_options system thread unit_test_framework
)
find_package(PETSC REQUIRED)
find_package(SPDLOG REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(LAPACK REQUIRED)
# Packages held in this repo in ./external
set(NLOHMANN_JSON_SINGLE_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/external/nlohmann_json")
set(GTEST_MPI_SINGLE_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/external/gtest-mpi-listener")
set(DEBUG_ASSERT_SINGLE_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/external/debug_assert")


set(CHYPS_SOURCE "./chyps")

# Setup necessary files to embed git information
set(PRE_CONFIGURE_FILE "${CHYPS_SOURCE}/git.cpp.in")
set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/chyps/git.cpp")
include(cmake/Utils/git_watcher.cmake)
add_library(git STATIC ${POST_CONFIGURE_FILE})
target_include_directories(git PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(git check_git)

add_library(chyps STATIC)
set_property(TARGET chyps PROPERTY CXX_STANDARD 17)
target_include_directories(chyps PUBLIC "./")
target_include_directories(chyps PUBLIC "${NLOHMANN_JSON_SINGLE_HEADER}")
target_include_directories(chyps PUBLIC "${DEBUG_ASSERT_SINGLE_HEADER}")
target_link_libraries(chyps PUBLIC MPI::MPI_CXX)
target_link_libraries(chyps PUBLIC Hypre::Hypre)
target_link_libraries(chyps PUBLIC Metis::Metis)
target_link_libraries(chyps PUBLIC Mfem::Mfem)
target_link_libraries(chyps PUBLIC Precice::Precice)
target_link_libraries(chyps PUBLIC
  Boost::boost
  Boost::filesystem
  Boost::log
  Boost::log_setup
  Boost::program_options
  Boost::system
  Boost::thread
  Boost::unit_test_framework
)
target_link_libraries(chyps PUBLIC Petsc::Petsc)
target_link_libraries(chyps PUBLIC Spdlog::Spdlog)
target_link_libraries(chyps PUBLIC adios2::cxx11_mpi)
target_link_libraries(chyps PUBLIC ${LAPACK_LIBRARIES})
target_link_libraries(chyps PUBLIC git)

add_executable(chyps_heat)
target_link_libraries(chyps_heat PRIVATE chyps)
set_property(TARGET chyps_heat PROPERTY CXX_STANDARD 17)

add_executable(chyps_initializer)
target_link_libraries(chyps_initializer PRIVATE chyps)
set_property(TARGET chyps_initializer PROPERTY CXX_STANDARD 17)


add_executable(dummy_precice)
target_link_libraries(dummy_precice PRIVATE chyps)
set_property(TARGET dummy_precice PROPERTY CXX_STANDARD 17)

add_subdirectory("${CHYPS_SOURCE}")

# Directory for test files
# Will automatically download GoogleTest and compile it to
# use for the testing framework.
if("${BUILD_TESTING}" MATCHES "ON")
  include(CTest)
  enable_testing()
  add_subdirectory(./tests)
endif()


