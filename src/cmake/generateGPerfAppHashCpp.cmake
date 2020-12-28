#
# Replace 'evaluate::[whatever]' by NULL
#
# * evaluate.h (class evaluate): New class implementing the actual evaluation
#   function (applications) for the equations in Qucs.
#
if (UNIX)
 add_custom_command(
   OUTPUT gperfappgen.h
   COMMAND
     ${SED_TOOL} -e 's/evaluate::[a-zA-Z0-9_]*/NULL/g' <
     ${CMAKE_CURRENT_SOURCE_DIR}/applications.h >
     ${CMAKE_CURRENT_BINARY_DIR}/gperfappgen.h
   DEPENDS ${applications.h})
else()
add_custom_command(
  OUTPUT gperfappgen.h
  COMMAND
    ${SED_TOOL} -e "s/evaluate::[a-zA-Z0-9_]*/NULL/g" <
    ${CMAKE_CURRENT_SOURCE_DIR}/applications.h >
    ${CMAKE_CURRENT_BINARY_DIR}/gperfappgen.h
  DEPENDS ${applications.h})
endif()

#
# Compile gperfappgen * used to generate gperf input file (used in qucsator)
#
set(gperf_SRC gperfappgen.cpp gperfappgen.h)
add_executable(gperfappgen ${gperf_SRC})
target_include_directories(gperfappgen PRIVATE ${qucsator_BINARY_DIR}
                            ${qucsator_app_BINARY_DIR}
                            ${qucsator_app_SOURCE_DIR}
			    ${qucsator_app_SOURCE_DIR}/math
			    )

#
# Run gperfappgen, pipe to gperf input to gperfapphash.gph
#
add_custom_command(
  OUTPUT gperfapphash.gph
  COMMAND gperfappgen > ${CMAKE_CURRENT_BINARY_DIR}/gperfapphash.gph
  DEPENDS gperfappgen
  COMMENT "gperfappgen gets executed"
  VERBATIM)

#
# Run gperf, create hash table. * -I, Include the necessary system include files
# at the beginning of the code. * -m, Perform multiple iterations to minimize
# generated table. * Replace '{""}' by '{"",0}; (why?)
#
if (UNIX)
 add_custom_command(
   OUTPUT gperfapphash.cpp
   COMMAND ${GPERF_TOOL} -I -m 8 ${CMAKE_CURRENT_BINARY_DIR}/gperfapphash.gph >
           temp.gperf
   COMMAND ${SED_TOOL} -e 's/{""},/{"",0},/g' < temp.gperf > ${CMAKE_CURRENT_BINARY_DIR}/gperfapphash.cpp
   DEPENDS gperfapphash.gph)
else()
add_custom_command(
  OUTPUT gperfapphash.cpp
  COMMAND ${GPERF_TOOL} -I -m 8 ${CMAKE_CURRENT_BINARY_DIR}/gperfapphash.gph >
          temp.gperf
  COMMAND ${SED_TOOL} -e "s/{\"\"},/{\"\",0},/g" < temp.gperf > ${CMAKE_CURRENT_BINARY_DIR}/gperfapphash.cpp
  DEPENDS gperfapphash.gph)
endif()
