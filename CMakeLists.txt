cmake_minimum_required(VERSION 3.7...3.18)

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


set(CHYPS_SOURCE "./chypss")

add_executable(chyps_heat)
set_property(TARGET chyps_heat PROPERTY CXX_STANDARD 11)
target_include_directories(chyps_heat PUBLIC "./")
target_link_libraries(chyps_heat PUBLIC MPI::MPI_CXX)
target_link_libraries(chyps_heat PUBLIC Hypre::Hypre)
target_link_libraries(chyps_heat PUBLIC Metis::Metis)
target_link_libraries(chyps_heat PUBLIC Mfem::Mfem)
target_link_libraries(chyps_heat PUBLIC Precice::Precice)
target_link_libraries(chyps_heat PUBLIC
  Boost::boost
  Boost::filesystem
  Boost::log
  Boost::log_setup
  Boost::program_options
  Boost::system
  Boost::thread
  Boost::unit_test_framework
)
target_link_libraries(chyps_heat PUBLIC Petsc::Petsc)
add_subdirectory("${CHYPS_SOURCE}")
