cmake_minimum_required(VERSION 3.1)

project(students_research)

enable_testing()

# Set global CMake properties
# -----------------------------------------------------------------------
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake targets")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(PROJECT_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release")
endif()

# Specify Include Directories
# -----------------------------------------------------------------------
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/doc_tracker
  ${CMAKE_CURRENT_SOURCE_DIR}/otsu_cq
)

# Set up default settings.
# First, set up default compiler options.
include(compiler)

# Include some helpers.
include(macro_helpers)

# Setting output directories
# -----------------------------------------------------------------------
se_set_generic_output_dirs(${CMAKE_SOURCE_DIR} ${BUILD_SUFFIX})

# Adding subdirectories
# -----------------------------------------------------------------------
add_subdirectory(doc_tracker)
add_subdirectory(otsu_cq)
