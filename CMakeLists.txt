# ~~~

#
# You can change the install location by running cmake like this: mkdir build;
# cd build cmake .. -DCMAKE_INSTALL_PREFIX=/new/install/prefix

# Notes:
#
# * By default, the prefix is "/usr/local"
# * Use -DADMSXML_DIR=[path] to give the path containing admsXml submodules (ADMS)
# * Bison can be pointed to with -DBISON_DIR=[path]

# Running with QtCreator Open this CMakeLists.txt with 'Open File or Project..."
# Provide CMake arguments: * Set admsXml path, (after building with Autotools) *
# Enable debug symbols for the debugger.
#
# Add the following arguments to the 'Run CMake': -DAMSXML_DIR=~/git/qucs/qucs-
# core/adms/admsXml/ -DCMAKE_BUILD_TYPE:STRING=Debug
#
# With Homebrew/MacPorts on Mac OS X: It might be necessary to change the PATH
# variable on the Build Environment. Prepend '/usr/local/bin' before '/usr/bin/'
# otherwise it will find the older executables provided by Apple.

# Build
# * qucsator
# * libqucs
# * qucsconv
#

# ~~~

project(qucs-core CXX)
cmake_minimum_required(VERSION 2.6)

# ignore the project() managed VERSION (new in CMake 3.0)
if(POLICY CMP0048)
  cmake_policy(SET CMP0048 OLD)
endif(POLICY CMP0048)


# otherwise qucsator cannot generate qucsdefs.h
add_definitions(-DDEBUG)

# defines nr_double_t
add_definitions(-DHAVE_CONFIG_H)

# ~~~
# TODO
# OPTION(ENABLE_QUCSLIB "enable qucslib build, default: OFF")
# OPTION(ENABLE_DOUBLE " type of double representation, default=double")

# TODO configure debug/release flags

# TODO check flags used on Autotools are needed:
#
# -pipe : Use pipes rather than temporary files for communication between the
# various stages of compilation.
#
# https://blog.mozilla.org/nnethercote/2011/01/18/the-dangers-of-fno-exceptions/
# -fno-exceptions option is used, which means that exception-handling is
# disabled.
#
# -fno-rtti :  Disable generation of information about every class with virtual
# functions.
#
# -fno-check-new : specific to GCC
#
# -Wmissing-prototypes : Warn if a global function is defined without a previous
# prototype declaration.
#
# ~~~

# use top VERSION file
file(STRINGS ${qucs-core_SOURCE_DIR}/VERSION QUCS_VERSION)
message(STATUS "Configuring ${PROJECT_NAME}: VERSION ${QUCS_VERSION}")

set(PROJECT_VERSION "${QUCS_VERSION}")

set(PROJECT_VENDOR "Qucs team. This program is licensed under the GNU GPL")

set(QUCS_URL "https://sf.net/p/qucs")
set(QUCS_BUGREPORT "qucs-bugs@lists.sourceforge.net")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/qucsconv.1.cmake.in
               ${CMAKE_CURRENT_SOURCE_DIR}/doc/qucsconv.1)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/qucsator.1.cmake.in
               ${CMAKE_CURRENT_SOURCE_DIR}/doc/qucsator.1)

# use last git commit hash along the version
set(GIT unknown)
if(EXISTS ${CMAKE_SOURCE_DIR}/../.git)
  find_package(Git)
  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log --pretty=format:%h -n 1u
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH)
  set(GIT ${GIT_COMMIT_HASH})
  message(STATUS "Found Git repository, last commit hash: ${GIT}")
endif()

# TODO rename the above variables? Project/Package?

# Define to the address where bug reports for this package should be sent.
set(PACKAGE_BUGREPORT "qucs-bugs@lists.sourceforge.net")

# Define to the full name of this package.
set(PACKAGE_NAME "qucs-core")

# Define to the full name and version of this package.
set(PACKAGE_STRING "${PACKAGE_NAME} ${PROJECT_VERSION}")

# Define to the one symbol short name of this package.
set(PACKAGE_TARNAME ${PACKAGE_NAME})

# Define to the home page for this package.
set(PACKAGE_URL "http://sourceforge.net/projects/qucs/")

# Define to the version of this package.
set(PACKAGE_VERSION ${PROJECT_VERSION})

#
# Avoid source tree pollution
#
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(
    FATAL_ERROR
      "\nIn-source builds are not permitted.
                       Make a separate folder for building:
                       $ mkdir build; cd build; cmake ..
                       Before that, remove the files already created:
                       $ rm -rf CMakeCache.txt CMakeFiles")
endif()

#
# Set locations of CMake modules, used on tests, find,...
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

#
# Need Flex
#
find_package(FLEX 2.5.9 REQUIRED)
if(FLEX_FOUND)
endif()

# ~~~
# Need Bison
#
# This is a HACK to get around a PATH issue with Qt Creator on OSX.
# It seams impossible to pass a custom PATH to Qt Creator on OSX, ie, cannot prepend `/usr/local/bin/` for instance.
# The FIND_PACKAGE fails. For now we provide a fallback with a custom FIND_PROGRAM.
# The variable BISON_DIR is also available.
# ~~~
if(WIN32)
  find_package(BISON 2.4 REQUIRED)
else() # Linux, OSX
  # use -DBISON_DIR=/path/ to provide the path to bison
  find_program(
    BISON_EXECUTABLE bison
    PATHS /usr/local/bin/ /opt/local/bin/ /usr/bin ${BISON_DIR}
    DOC "bison path"
    NO_DEFAULT_PATH)
  if(BISON_EXECUTABLE)
    message(STATUS "Found bison: " ${BISON_EXECUTABLE})
  else()
    message(
      FATAL_ERROR "Unable to find bison. Try to provide -DBISON_DIR=[path]")
  endif()
endif()

#
# Check for sed
#
find_program(SED_TOOL NAMES sed)
if(NOT SED_TOOL)
  message(FATAL_ERROR "Unable to find sed")
else()
  message(STATUS "Found sed: " ${SED_TOOL})
endif()

#
# Check for gperf
#
find_program(GPERF_TOOL NAMES gperf)
if(NOT GPERF_TOOL)
  message(FATAL_ERROR "gperf required in PATH")
else()
  message(STATUS "Found gperf: " ${GPERF_TOOL})
endif()

#
# Check if admsXml is available
#
# * Use -DADMSXML_DIR=[path] to give the path containing admsXml
# * Try a few other locations
#
find_program(
  ADMSXML admsXml
  HINTS ${ADMSXML_DIR}
  PATHS /usr/local/bin/ /opt/local/bin/ /usr/bin
  DOC "admsXml application")
if(NOT ADMSXML)
  message(FATAL_ERROR "admsXml required in PATH")
else()
  message(STATUS "Found admsXml: " ${ADMSXML})
endif()

#
# Set up RPATH for the project
#
option(ENABLE_RPATH "Enable rpath support on Linux and Mac" ON)
if(NOT CMAKE_INSTALL_RPATH)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()
if(APPLE AND NOT CMAKE_INSTALL_NAME_DIR)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
endif()
if(UNIX AND ENABLE_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_MACOSX_RPATH TRUE)
endif()

# CMake adds --enable-all-exports on Cygwin (since Cygwin is supposed to be
# UNIX-like), but we need to add it explicitly for a native windows build with
# the MinGW tools.
if(WIN32)
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
      "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
      "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-auto-import")
endif()

# ~~~
# indiscriminate copy/paste from:
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake/20165220#20165220

# Initialize CXXFLAGS.

# \todo fix headers and use standard C++ methods * strdup is not C or C++
# standard, it is POSIX adding -fpermissive let it compile with a ton of
# warnings * problem with non-starndart _stricmp using -stdr=c++0x set g++ into
# strict ANSY, relax that with -U__STRICT_ANSI__. Could use -std=gnu++0x
# ~~~
if(WIN32)
  set(CMAKE_CXX_FLAGS "-Wall -std=c++0x -fpermissive -U__STRICT_ANSI__")
else()
  set(CMAKE_CXX_FLAGS "-Wall -std=c++11")
endif()

# indiscriminate copy/paste from:
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-
# compiler-with-cmake/20165220#20165220
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# Compiler-specific C++11 activation.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  if(NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
    message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
  endif()
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else()
  message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif()

#
# Set position independent code PIC
#
if(UNIX AND NOT APPLE)
  if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  endif()
endif()

#
# Go look for stuff to build/install...
#
add_subdirectory(src)
add_subdirectory(doc)

#
# Custom uninstall target
#
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake" IMMEDIATE @ONLY)

add_custom_target(
  uninstall-core COMMAND ${CMAKE_COMMAND} -P
                        ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake)

# TODO install distributables EXTRA_DIST = BUGS bootstrap depcomp RELEASE

# TODO tarball TODO bundle

set(CPACK_GENERATOR "TGZ")
# SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "guitorri") #required

# build a CPack driven installer package include
# (InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_VERSION_MAJOR 0)
set(CPACK_PACKAGE_VERSION_MINOR 18)
include(CPack)
