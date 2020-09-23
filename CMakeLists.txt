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


add_executable(heat)
set_property(TARGET heat PROPERTY CXX_STANDARD 11)
target_include_directories(heat PRIVATE /usr/local/include)
target_include_directories(heat PRIVATE /Users/rmc298/lib_builds/hypre/include)
target_include_directories(heat PRIVATE "${MPI_INCLUDE_PATH}")
target_link_libraries(heat PRIVATE "${MPI_CXX_LIBRARIES}")
target_link_libraries(heat PRIVATE /usr/local/lib/libmfem.a)
target_link_libraries(heat PRIVATE /Users/rmc298/lib_builds/hypre/lib/libHYPRE.a)
target_link_libraries(heat PRIVATE /Users/rmc298/lib_builds/hypre/lib/libHYPRE.a)
target_link_libraries(heat PRIVATE /usr/local/lib/libmetis.dylib)

add_executable(euler)
set_property(TARGET euler PROPERTY CXX_STANDARD 11)
target_include_directories(euler PRIVATE /usr/local/include)
target_include_directories(euler PRIVATE /Users/rmc298/lib_builds/hypre/include)
target_include_directories(euler PRIVATE "${MPI_INCLUDE_PATH}")
target_link_libraries(euler PRIVATE "${MPI_CXX_LIBRARIES}")
target_link_libraries(euler PRIVATE /usr/local/lib/libmfem.a)
target_link_libraries(euler PRIVATE /Users/rmc298/lib_builds/hypre/lib/libHYPRE.a)
target_link_libraries(euler PRIVATE /Users/rmc298/lib_builds/hypre/lib/libHYPRE.a)
target_link_libraries(euler PRIVATE /usr/local/lib/libmetis.dylib)

add_subdirectory("${CHYPSS_SOURCE}")
