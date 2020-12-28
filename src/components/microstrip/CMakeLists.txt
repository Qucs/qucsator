project(microstrip)
cmake_minimum_required(VERSION 3.16)

include(${PROJECT_SOURCE_DIR}/../../../../cmake/readMakeList.txt)

#qucsator_SOURCE_DIR is defined as soon as a project named qucsator is created
if ((NOT qucsator_BINARY_DIR) OR (NOT qucsator_app_SOURCE_DIR) )
	message(FATAL_ERROR "This project cannot be used as standalone, because qucs_typedef is missing. Everything around the command
		configure_file(\"${qucsator_SOURCE_DIR}/qucs_typedefs.h.cmake\"
			\"${qucsator_BINARY_DIR}/qucs_typedefs.h\")\
		is needed to get it work. But this must be done in all subfiles and this makes no sense.")
endif()

readVariable(MakeList SRCS MICROSTRIP_SRC)
readVariable(MakeList HDRS MICROSTRIP_HDRS)

add_library(coreMicrostrip OBJECT ${MICROSTRIP_SRC})
target_include_directories(coreMicrostrip
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} # This folder
	PRIVATE
	${qucsator_BINARY_DIR} # qucs_typedefs.h and config.h
	${qucsator_app_SOURCE_DIR}
	${qucsator_app_SOURCE_DIR}/components
	${qucsator_app_SOURCE_DIR}/components/devices
	${qucsator_app_SOURCE_DIR}/math)

install(FILES ${MICROSTRIP_HDRS} DESTINATION include/qucsator)

