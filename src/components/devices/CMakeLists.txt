project(devices)
cmake_minimum_required(VERSION 3.16)

include(${PROJECT_SOURCE_DIR}/../../../../cmake/readMakeList.txt)

#qucsator_SOURCE_DIR is defined as soon as a project named qucsator is created
if ((NOT qucsator_BINARY_DIR) OR (NOT qucsator_app_SOURCE_DIR) )
	message(FATAL_ERROR "This project cannot be used as standalone, because qucs_typedef is missing. Everything around the command
		configure_file(\"${qucsator_SOURCE_DIR}/qucs_typedefs.h.cmake\"
			\"${qucsator_BINARY_DIR}/qucs_typedefs.h\")\
		is needed to get it work. But this must be done in all subfiles and this makes no sense.")
endif()

readVariable(MakeList SRCS DEVICES_SRC)
readVariable(MakeList HDRS DEVICES_HEADERS)

add_library(coreDevices OBJECT ${DEVICES_SRC})
target_include_directories(coreDevices
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} # This folder
	PRIVATE
	${qucsator_BINARY_DIR}
	${qucsator_app_SOURCE_DIR}
	${qucsator_app_SOURCE_DIR}/components
	${qucsator_app_SOURCE_DIR}/math)

install(FILES ${DEVICES_HEADERS} DESTINATION include/qucsator)
