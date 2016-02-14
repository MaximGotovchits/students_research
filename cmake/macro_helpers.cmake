# Here we put various helpers for CMakeLists.txt scripts.

INCLUDE(target_arch)



MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
        SET(dirlist ${dirlist} ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()



FUNCTION(SE_DETECT_PLATFORM platform)
  STRING(TOLOWER ${CMAKE_SYSTEM_NAME} system)
  IF (${system} MATCHES "generic" AND CMAKE_CROSSCOMPILING)
    SET(system "${TOOLCHAIN_SYSTEM_NAME}")
  ENDIF()

  target_architecture(architecture)

  SET(${platform} "${architecture}-${system}" PARENT_SCOPE)
ENDFUNCTION()


# SE_SET_GENERIC_OUTPUT_DIRS( root_dir [build_suffix] )
MACRO(SE_SET_GENERIC_OUTPUT_DIRS root_dir)
  # Determine if optional BUILD_SUFFIX is presented

  # Cannot use ARGN directly with list() command.
  # Copy to a variable first.
  SET(__extra_macro_args ${ARGN})
  # Did we get any optional args?
  LIST(LENGTH __extra_macro_args __num_extra_args)
  IF (${__num_extra_args} GREATER 0)
    LIST(GET __extra_macro_args 0 __build_suffix)
  ENDIF ()

  SE_DETECT_PLATFORM(__platform)
  message(STATUS "Platform is ${__platform}")

  # -----------------------------------------------------------------------
  # Runtime output directory
  # -----------------------------------------------------------------------
  SET(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${root_dir}/bin.${__platform}.${CMAKE_BUILD_TYPE}")
  SET(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG
    "${root_dir}/bin.${__platform}.debug")
  SET(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE
    "${root_dir}/bin.${__platform}.release")
  SET(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO
    "${root_dir}/bin.${__platform}.relwithdebinfo")

  IF (__build_suffix)
    SET(
      CMAKE_RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}.${__build_suffix}")
    SET(
      CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}.${__build_suffix}")
    SET(
      CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}.${__build_suffix}")
    SET(
      CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO}.${__build_suffix}")
  ENDIF()

  # -----------------------------------------------------------------------
  # Set library output directories
  # -----------------------------------------------------------------------
  SET(
    CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_SOURCE_DIR}/bin.${__platform}.${CMAKE_BUILD_TYPE}")
  SET(
    CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG
    "${CMAKE_SOURCE_DIR}/bin.${__platform}.debug")
  SET(
    CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE
    "${CMAKE_SOURCE_DIR}/bin.${__platform}.release")
  SET(
    CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO
    "${CMAKE_SOURCE_DIR}/bin.${__platform}.relwithdebinfo")
  
  IF (__build_suffix)
    SET(
      CMAKE_LIBRARY_OUTPUT_DIRECTORY
      "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}.${__build_suffix}")
    SET(
      CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG
      "${CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG}.${__build_suffix}")
    SET(
      CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE
      "${CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE}.${__build_suffix}")
    SET(
      CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO
      "${CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO}.${__build_suffix}")
  ENDIF()

  # -----------------------------------------------------------------------
  # Set archive output directory
  # -----------------------------------------------------------------------
  SET(
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_SOURCE_DIR}/lib/${__platform}.${CMAKE_BUILD_TYPE}")
  SET(
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG
    "${CMAKE_SOURCE_DIR}/lib/${__platform}.debug")
  SET(
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE
    "${CMAKE_SOURCE_DIR}/lib/${__platform}.release")
  SET(
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO
    "${CMAKE_SOURCE_DIR}/lib/${__platform}.relwithdebinfo")
  IF (__build_suffix)
    SET(
      CMAKE_ARCHIVE_OUTPUT_DIRECTORY
      "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}.${__build_suffix}")
    SET(
      CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG
      "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG}.${__build_suffix}")
    SET(
      CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE
      "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE}.${__build_suffix}")
    SET(
      CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO
      "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO}.${__build_suffix}")
  ENDIF()

  # -----------------------------------------------------------------------
  # Set install directories
  # -----------------------------------------------------------------------
  SET(CMAKE_INSTALL_PREFIX "${PROJECT_ROOT_DIR}/install.${__platform}.${CMAKE_BUILD_TYPE}")
ENDMACRO()