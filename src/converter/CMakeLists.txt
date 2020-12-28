project(converter)
cmake_minimum_required(VERSION 3.16)

include(${PROJECT_SOURCE_DIR}/../../../cmake/readMakeList.txt)

find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)

#include_directories(
#  ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/math
#  ${CMAKE_CURRENT_BINARY_DIR}) # qucdefs.h

readVariable(MakeList SRCS QUCSCONV_SRC)

set(ParserTypes spice vcd)

set(conv_generated)
foreach(type ${ParserTypes})
  # Create custom Bison
  set(bisonIn "${CMAKE_CURRENT_SOURCE_DIR}/parse_${type}.ypp")
  set(bisonOut "parse_${type}.hpp"
               "parse_${type}.cpp")
  add_custom_command(
    OUTPUT ${bisonOut}
    COMMAND
      ${BISON_EXECUTABLE}
      --defines=parse_${type}.hpp
      --output=parse_${type}.cpp
      ${bisonIn}
    DEPENDS ${bisonIn})
  # Create custom Flex
  set(flexIn "${CMAKE_CURRENT_SOURCE_DIR}/scan_${type}.lpp")
  set(flexOut "scan_${type}.cpp")
  add_custom_command(
    OUTPUT ${flexOut}
    COMMAND ${FLEX_EXECUTABLE} --outfile=${flexOut} ${flexIn}
    DEPENDS ${flexIn})
  list(APPEND conv_generated ${bisonOut})
  list(APPEND conv_generated ${flexOut})
endforeach()

#
# qucsator is needed before qucsconv can be compiled run qucsator -l > pipe to
# qucsdef.h
#
add_custom_command(
  OUTPUT qucsdefs.h
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../qucsator -l >
          ${CMAKE_CURRENT_BINARY_DIR}/qucsdefs.h
  DEPENDS qucsator) # TODO: here is a build problem

add_custom_target(defs_converter DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/qucsdefs.h)

# TODO: check that the number of files in conv_generated is same as for ParserTypes
add_executable(qucsconv ${QUCSCONV_SRC} ${conv_generated})

add_dependencies(qucsconv defs_converter)

# TODO: why CMAKE_DL_LIBS needed?
target_link_libraries(qucsconv libqucsator ${CMAKE_DL_LIBS})
target_include_directories(qucsconv
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
	PRIVATE ${qucsator_BINARY_DIR})

#
# Handle installation
#
install(TARGETS qucsconv DESTINATION bin)
