project(interface)
cmake_minimum_required(VERSION 3.16)

include(${PROJECT_SOURCE_DIR}/../../../cmake/readMakeList.txt)


#include_directories(
#  ${qucs-core_SOURCE_DIR}
#  ${qucs-core_CURRENT_SOURCE_DIR}
#  ${qucs-core_SOURCE_DIR}/src/math
#  ${qucs-core_SOURCE_DIR}/src/components # component.h
#  ${qucs-core_BINARY_DIR}/src/components # generated verilog/[].core.h
#)

readVariable(MakeList SRCS INTERFACE_SRC)
readVariable(MakeList HDRS HEADERS)

add_library(coreInterface OBJECT ${INTERFACE_SRC})
target_include_directories(coreInterface
	PUBLIC
	${qucsator_app_BINARY_DIR}/components  # verilog/[].core.h
	${CMAKE_CURRENT_SOURCE_DIR}
	PRIVATE
	${qucsator_BINARY_DIR} # config.h
	${qucsator_app_SOURCE_DIR} # compat.h
	${qucsator_app_SOURCE_DIR}/math # precision.h
	${qucsator_app_SOURCE_DIR}/components # component.h
	)

# coreVerilog is located in components/verilog
add_dependencies(coreInterface coreVerilog)

install(FILES ${HEADERS} DESTINATION include/qucsator)
