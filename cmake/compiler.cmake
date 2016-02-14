# Default compiler flags.

# Flags for MS Visual Studio.
IF (MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_DEPRECATE)

  # Classic workaround to prevent windows.h replace std::min/std::max with
  # its macros.
  ADD_DEFINITIONS(-DNOMINMAX)

  # Workaround MSVC bug with template specialization storage type.
  ADD_DEFINITIONS(-DUSE_STATIC_SPECIAL)

# It is GCC or Clang.
ELSEIF (CMAKE_COMPILER_IS_GNUCXX OR
        CMAKE_COMPILER_IS_GNUC OR
        ${CMAKE_C_COMPILER_ID} MATCHES "Clang" OR
        ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

  # Make all libc functions re-entrant.
  ADD_DEFINITIONS(-D_REENTRANT)

  # Hide emscripten warnings
  IF (EMSCRIPTEN)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-warn-absolute-paths")
  ENDIF()

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wcast-align")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused -Wno-missing-field-initializers")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-cast-align")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-attributes")
  IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUC)
     # See http://stackoverflow.com/questions/6687630/c-c-gcc-ld-remove-unused-symbols
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffunction-sections -fdata-sections")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--gc-sections")
     set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--gc-sections")
  ELSEIF(${CMAKE_C_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
     # http://stackoverflow.com/questions/17710024/clang-removing-dead-code-during-static-linking-gcc-equivalent-of-wl-gc-sect
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -dead_strip")
     set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -dead_strip")
  ENDIF()

  IF (${CMAKE_SYSTEM_PROCESSOR} MATCHES "^arm.*$")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfloat-abi=softfp")
  ENDIF()

  # Debug
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS}")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -ggdb -O0 -g3 -fno-inline -fno-omit-frame-pointer")
  # Release
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS}")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Os")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -DNDEBUG")
  # Release with debug info (for profile purposes).
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${CMAKE_C_FLAGS_RELEASE}")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g -fno-omit-frame-pointer")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG")

  # CXX stuff is same as C stuff, plus something C++ specific.
  SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fexceptions")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fexceptions")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fexceptions")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -fexceptions")

  # enable C++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# It is unknown. Do nothing.
ELSE()

  MESSAGE(WARNING "Unknown compiler. Not using any specific flags. The build will probably fail.")

ENDIF()
