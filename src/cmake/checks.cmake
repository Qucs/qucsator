#
# Checks for libraries.
#
# AC_CHECK_LIB(m, sin) need to check for sin?
if(NOT WIN32)
  find_library(MATH_LIB NAMES m)
  if(NOT MATH_LIB)
	  message(SEND_ERROR "Math lib not found: ${MATH_LIB}")
else()
	  message(STATUS "Math lib found at: ${MATH_LIB}")
endif()
endif()

#
# Checks for header files. AC_HEADER_STDC !! obsolete, need to check? Define
# STDC_HEADERS if the system has ANSI C header files. Specifically, this macro
# checks for `stdlib.h', `stdarg.h', `string.h', and `float.h'; Lifted the cmake
# checks from gd-libdg, Lua https://bitbucket.org/libgd/gd-libgd
# https://github.com/LuaDist/libgd/tree/master/cmake/modules
include(CheckIncludeFiles)
set(CMAKE_REQUIRED_INCLUDES "/usr/include" "/usr/local/include")
set(CMAKE_MODULE_PATH "${qucsator_SOURCE_DIR}/cmake/modules")
include(AC_HEADER_STDC)

#
# Further header checks
#
include(CheckIncludeFile)

# list of headers to be checked
set(INCLUDES ieeefp.h memory.h stddef.h stdlib.h string.h unistd.h)

#
# Check if header can be included. * Define HAVE_[basename]_H to 1 if you have
# the header.
#
foreach(header ${INCLUDES})
  get_filename_component(base ${header} NAME_WE)
  string(TOUPPER ${base} base)
  check_include_file(${header} HAVE_${base}_H)
  # MESSAGE(STATUS "${header}  --> ${HAVE_${base}_H}")
endforeach()

# Checks for typedefs, structures, and compiler characteristics. AC_C_CONST
# !!obsolete AC_C_CONST "This macro is obsolescent, as current C compilers
# support `const'. New programs need not use this macro."

#
# Check for library functions * not all functons seem to be used after defined.
# TODO check for HAVE_{func}
#
include(CheckFunctionExists)
set(REQUIRED_FUNCTIONS
	floor
	pow
	exp
	sqrt
	log10
	log
	cos
	sin
	acos
	asin # for real.cpp
	tan
	atan
	sinh
	cosh
	tanh
	fabs
	modf
	atan2
	jn
	yn
	erf
	erfc # for fspecial.cpp
	round
	trunc
	acosh
	asinh # for real.cpp
	strdup
	strerror
	strchr) # for compat.h, matvec.cpp, scan_*.cpp

foreach(func ${REQUIRED_FUNCTIONS})
  string(TOUPPER ${func} FNAME)
  check_function_exists(${func} HAVE_${FNAME})
  # message(STATUS "${func}  --> ${HAVE_${FNAME}}")
endforeach()

#
# Checks for complex classes and functions, as in the Autotools scripts.
#
# AC_CXX_NAMESPACES !!custom m4 AC_CXX_HAVE_COMPLEX !!custom m4
# AC_CXX_HAVE_TR1_COMPLEX !!custom m4 AC_CHECK_CXX_COMPLEX_FUNCS([cos cosh exp
# log log10 sin sinh sqrt tan tanh]) !!custom m4
# AC_CHECK_CXX_COMPLEX_FUNCS([acos acosh asin asinh atan atanh])
# AC_CHECK_CXX_COMPLEX_FUNCS([log2 norm]) AC_CHECK_CXX_COMPLEX_POW
# AC_CHECK_CXX_COMPLEX_ATAN2         !failed, need libstdc?
# AC_CHECK_CXX_COMPLEX_FMOD          !failed: AC_CHECK_CXX_COMPLEX_POLAR
# AC_CHECK_CXX_COMPLEX_POLAR_COMPLEX !failed

#
# Namespace
#
# Check whether the compiler implements namespaces
#
try_compile(
  HAVE_NAMESPACES ${CMAKE_BINARY_DIR}
  ${qucsator_SOURCE_DIR}/cmake/namespaces.cpp OUTPUT_VARIABLE TRY_OUT)
if(NOT HAVE_NAMESPACES)
  message(
	  SEND_ERROR
	"${PROJECT_NAME} requires an c++ compiler with namespace HAVE_NAMESPACES failed"
  ) # ${TRY_OUT}")
endif()

#
# Check whether the compiler has complex<T>
#
try_compile(HAVE_COMPLEX ${CMAKE_BINARY_DIR}
	${qucsator_SOURCE_DIR}/cmake/complex.cpp OUTPUT_VARIABLE TRY_OUT)
if(NOT HAVE_COMPLEX)
  message(SEND_ERROR "HAVE_COMPLEX failed") # ${TRY_OUT}")
endif()

#
# Check std::vector::erase iterator type [const_iterator | iterator]
#
message(STATUS "Checking HAVE_ERASE_CONSTANT_ITERATOR")
try_compile(
  HAVE_ERASE_CONSTANT_ITERATOR ${CMAKE_BINARY_DIR}
  ${qucsator_SOURCE_DIR}/cmake/erase_iterator_type.cpp OUTPUT_VARIABLE TRY_OUT)
if(HAVE_ERASE_CONSTANT_ITERATOR)
  message(STATUS "Using std::vector:erase iterator type : const_iterator")
else()
  message(STATUS "Using std::vector:erase iterator type : iterator")
endif()

#
# Check for list of complex functions.
#
set(COMPLEX_FUNCS_GRP1
	acos
	acosh
	asin
	asinh
	atan
	atanh
	cos
	cosh
	exp
	log
	log10
	sin
	sinh
	sqrt
	tan
	tanh
	log2
	norm)
set(COMPLEX_FUNCS_GRP2 pow atan2 fmod polar polar_complex)

#
# test complex function group 1 code inlined to easily replace '${func}'
#
foreach(func ${COMPLEX_FUNCS_GRP1})
  set(code
	  " #include <complex>
	  using namespace std\;
	  #ifdef log2
	#undef log2
	#endif

	int main() {
	complex<double> a\;
	  ${func}(a)\;
	  return 0\;
	  }")

file(WRITE ${qucsator_SOURCE_DIR}/cmake/test_${func}.cpp ${code})

  string(TOUPPER ${func} FNAME)

  message(STATUS "Checking HAVE_CXX_COMPLEX_${FNAME}")

  try_compile(
	  HAVE_CXX_COMPLEX_${FNAME} ${CMAKE_BINARY_DIR}
	${qucsator_SOURCE_DIR}/cmake/test_${func}.cpp OUTPUT_VARIABLE TRY_OUT)
if(NOT HAVE_CXX_COMPLEX_${FNAME})
	  message(STATUS "HAVE_CXX_COMPLEX_${FNAME} failed") # ${TRY_OUT}")
endif()

  file(REMOVE ${qucsator_SOURCE_DIR}/cmake/test_${func}.cpp ${code})
endforeach()

#
# test complex function group 2 use prepared source file.
#
foreach(func ${COMPLEX_FUNCS_GRP2})
  string(TOUPPER ${func} FNAME)

  message(STATUS "Checking HAVE_CXX_COMPLEX_${FNAME}")

  try_compile(
	  HAVE_CXX_COMPLEX_${FNAME} ${CMAKE_BINARY_DIR}
	${qucsator_SOURCE_DIR}/cmake/complex_${func}.cpp
	COMPILE_DEFINITIONS -DHAVE_NAMESPACES -DHAVE_COMPLEX -DHAVE_TR1_COMPLEX
	OUTPUT_VARIABLE TRY_OUT)
if(NOT HAVE_CXX_COMPLEX_${FNAME})
	  message(STATUS "HAVE_CXX_COMPLEX_${FNAME} failed") # ${TRY_OUT}")
endif()
endforeach()
