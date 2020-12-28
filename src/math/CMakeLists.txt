project(math)
cmake_minimum_required(VERSION 3.16)

include(${PROJECT_SOURCE_DIR}/../../../cmake/readMakeList.txt)

#qucsator_SOURCE_DIR is defined as soon as a project named qucsator is created
if ((NOT qucsator_BINARY_DIR) OR (NOT qucsator_app_SOURCE_DIR) )
	message(FATAL_ERROR "This project cannot be used as standalone, because qucs_typedef is missing. Everything around the command
		configure_file(\"${qucsator_SOURCE_DIR}/qucs_typedefs.h.cmake\"
			\"${qucsator_BINARY_DIR}/qucs_typedefs.h\")\
		is needed to get it work. But this must be done in all subfiles and this makes no sense.")
endif()

readVariable(MakeList SRCS MATH_SRC)
readVariable(MakeList HDRS HEADERS)

add_library(coreMath OBJECT ${MATH_SRC})
# ${qucsator_BINARY_DIR} must be added, because from the qucsator/qucs_typedefs.h.cmake the header is created and this one is stored inside the Binary directory
# constants.h from qucs/src is needed in the math library. So the path must also be added
target_include_directories(coreMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(coreMath PRIVATE
                           ${qucsator_app_SOURCE_DIR}
                           ${qucsator_BINARY_DIR} # config.h and qucs_typedefs.h
)

install(FILES ${HEADERS} DESTINATION include/qucsator)
