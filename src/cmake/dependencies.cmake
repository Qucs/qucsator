# In this file all dependencies of libqucsator are listed. This file gets exectued from the libqucsator CMakeLists.txt

find_program(SED_TOOL NAMES sed REQUIRED)
message(STATUS "Found sed: " ${SED_TOOL})

find_program(GPERF_TOOL NAMES gperf REQUIRED)
message(STATUS "Found gperf: " ${GPERF_TOOL})

#
# Need Flex
#
find_package(FLEX 2.5.4 REQUIRED)

# ~~~
# Need Bison
#
# This is a HACK to get arround a PATH issue with Qt Creator on OSX.
# It seams impossible to pass a custom PATH to Qt Creator on OSX, ie, cannot prepend `/usr/local/bin/` for intance.
# The FIND_PACKAGE fails. For now we provide a fallback with a custom FIND_PROGRAM.
# The variable BISON_DIR is also available.
# ~~~
if(WIN32 OR (UNIX AND NOT APPLE))
  find_package(BISON 2.4 REQUIRED)
elseif(APPLE) # OSX. TODO: check if still relevant
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
# Check if admsXml is available
#
# * Use -DADMSXML_DIR=[path] to give the path containing admsXml
# * Try a few other locations
#
# Needed by src/components/verlog
# https://stackoverflow.com/questions/36102809/qucs-core-configure-error-needs-admsxml
# TODO: download if not already installed
find_program(
  ADMSXML admsXml
  HINTS ${ADMSXML_DIR}
  PATHS /usr/local/bin/ /opt/local/bin/ /usr/bin
  REQUIRED
  DOC "admsXml application")
message(STATUS "Found admsXml: " ${ADMSXML})
