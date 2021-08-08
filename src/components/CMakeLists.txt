project(components)
cmake_minimum_required(VERSION 3.16)

include(${PROJECT_SOURCE_DIR}/../../../cmake/readMakeList.txt)

readVariable(MakeList SRCS COMPONENTS_SRC)
readVariable(MakeList HDRS COMPONENTS_HDRS)

add_library(coreComponents OBJECT ${COMPONENTS_SRC})
target_include_directories(coreComponents
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} # This folder
	PRIVATE ${qucsator_BINARY_DIR} # config.h and qucs_typedef.h
        ${qucsator_app_SOURCE_DIR}
		${qucsator_app_SOURCE_DIR}/math
	)

install(FILES ${COMPONENTS_HDRS} DESTINATION include/qucsator)
