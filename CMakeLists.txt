cmake_minimum_required(VERSION 3.7...3.18)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(CHyPSS VERSION 0.1
               DESCRIPTION "MFEM Learning Code"
               LANGUAGES CXX)

find_package(MPI)
message("${MPI_CXX_LIBRARIES}")


set(CHYPSS_SOURCE "./chypss")
set(HYPRE_DIR "/Users/rmc298/lib_builds/hypre")
set(METIS_DIR "/usr/local/")
set(MFEM_DIR "/usr/local/")


add_executable(heat)
set_property(TARGET heat PROPERTY CXX_STANDARD 11)
target_include_directories(heat PRIVATE "${METIS_DIR}/include")
target_include_directories(heat PRIVATE "${MFEM_DIR}/include")
target_include_directories(heat PRIVATE "${HYPRE_DIR}/include")
target_include_directories(heat PRIVATE "${MPI_INCLUDE_PATH}")
target_link_libraries(heat PRIVATE "${MPI_CXX_LIBRARIES}")
target_link_libraries(heat PRIVATE "${MFEM_DIR}/lib/libmfem.a")
target_link_libraries(heat PRIVATE "${HYPRE_DIR}/lib/libHYPRE.a")
target_link_libraries(heat PRIVATE "${METIS_DIR}/lib/libmetis.dylib")

add_executable(euler)
set_property(TARGET euler PROPERTY CXX_STANDARD 11)
target_include_directories(euler PRIVATE "${METIS_DIR}/include")
target_include_directories(euler PRIVATE "${MFEM_DIR}/include")
target_include_directories(euler PRIVATE "${HYPRE_DIR}/include")
target_include_directories(euler PRIVATE "${MPI_INCLUDE_PATH}")
target_link_libraries(euler PRIVATE "${MPI_CXX_LIBRARIES}")
target_link_libraries(euler PRIVATE "${MFEM_DIR}/lib/libmfem.a")
target_link_libraries(euler PRIVATE "${HYPRE_DIR}/lib/libHYPRE.a")
target_link_libraries(euler PRIVATE "${METIS_DIR}/lib/libmetis.dylib")

add_subdirectory("${CHYPSS_SOURCE}")
