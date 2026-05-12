cmake_minimum_required(VERSION 3.20)
project(rigid_body_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(RBE_BUILD_TESTS "Build GoogleTest suite" ON)
option(RBE_BUILD_BENCH "Build benchmark CLI" ON)
option(RBE_BUILD_DEMO "Build the GLFW demo (OpenGL + Vulkan)" OFF)
option(RBE_TRACY "Instrument with Tracy zones" OFF)
option(RBE_SANITIZE "Build with AddressSanitizer and UBSan" OFF)

include(FetchContent)

add_library(rbe STATIC
  src/world.cpp
  src/collision.cpp
  src/dynamic_tree.cpp)
target_include_directories(rbe PUBLIC include PRIVATE src)

# Determinism: no fused multiply-add contraction, so the same build gives bit-identical results.
if(MSVC)
  target_compile_options(rbe PRIVATE /W4 /fp:precise)
else()
  target_compile_options(rbe PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
endif()

if(RBE_SANITIZE AND NOT MSVC)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
  target_compile_options(rbe PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(rbe PUBLIC -fsanitize=address,undefined)
endif()

if(RBE_TRACY)
  FetchContent_Declare(tracy GIT_REPOSITORY https://github.com/wolfpld/tracy.git GIT_TAG v0.11.1 GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(tracy)
  target_link_libraries(rbe PUBLIC Tracy::TracyClient)
  target_compile_definitions(rbe PUBLIC RBE_TRACY)
endif()

add_library(rbe_scenes STATIC scenes/scenes.cpp)
target_include_directories(rbe_scenes PUBLIC scenes)
target_link_libraries(rbe_scenes PUBLIC rbe)

if(RBE_BUILD_TESTS)
  enable_testing()
  find_package(GTest QUIET)
  if(NOT GTest_FOUND)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG v1.15.2 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(googletest)
  endif()
  add_executable(rbe_tests
    tests/test_math.cpp
    tests/test_collision.cpp
    tests/test_tree.cpp
    tests/test_physics.cpp)
  target_include_directories(rbe_tests PRIVATE src)
  target_link_libraries(rbe_tests PRIVATE rbe rbe_scenes GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(rbe_tests DISCOVERY_TIMEOUT 60)
endif()

if(RBE_BUILD_BENCH)
  add_executable(rbe_bench bench/rbe_bench.cpp)
  target_link_libraries(rbe_bench PRIVATE rbe rbe_scenes)
endif()

if(RBE_BUILD_DEMO)
  add_subdirectory(demo)
endif()
