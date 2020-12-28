# Generate Lexer and Parser files with Bison
#
# List of lexer/parsers type names
#
set(ParserTypes
        csv
	citi
	dataset
	mdl
	netlist
	touchstone
	zvr)

set(generated_SRC)
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
	    #COMMAND ${FLEX_EXECUTABLE} --outfile=${flexOut} ${flexIn} # not working on windows 10?
	    COMMAND ${FLEX_EXECUTABLE} -o"${flexOut}" ${flexIn}
	    DEPENDS ${flexIn})

    list(APPEND generated_SRC ${bisonOut})
    list(APPEND generated_SRC ${flexOut})
endforeach()
