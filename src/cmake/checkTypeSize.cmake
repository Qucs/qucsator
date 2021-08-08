#
# Check for type sizes.
#
include(CheckTypeSize)
check_type_size("short" SIZEOF_SHORT)
check_type_size("int" SIZEOF_INT)
check_type_size("long" SIZEOF_LONG)
check_type_size("double" SIZEOF_DOUBLE)
check_type_size("long double" SIZEOF_LONG_DOUBLE)
# MESSAGE(STATUS "short  ${SIZEOF_SHORT}" ) MESSAGE(STATUS "int   ${SIZEOF_INT}"
# ) MESSAGE(STATUS "long  ${SIZEOF_LONG}" ) MESSAGE(STATUS "double
# ${SIZEOF_DOUBLE}" ) MESSAGE(STATUS "long double ${SIZEOF_LONG_DOUBLE}" )

#
# Check for double type. * valid types are: double, float and long double. *
# defines: nr_double_t,  The global type of double representation. * defines:
# NR_DOUBLE_SIZE,  The size of the double representation. * Use -DENABLE-
# DOUBLE="[float,double, long double]"
if(ENABLE-DOUBLE)
  # User defined
  set(DoubleType ${ENABLE-DOUBLE})

  # valid types
  set(ValidTypes "float" "double" "long double")

  list(FIND ValidTypes ${DoubleType} HasType)
  if(HasType EQUAL -1)
          message(FATAL_ERROR "Valid types are: ${ValidTypes}")
  endif()

  # The global type of double representation.
  set(nr_double_t DoubleType)
  check_type_size(${DoubleType} DoubleSize)

  # The size of the double representation.
  set(NR_DOUBLE_SIZE ${DoubleSize})
else()
  # Default double
  set(DoubleType "double")
  # The global type of double representation.
  set(nr_double_t ${DoubleType})
  check_type_size(${DoubleType} DoubleSize)

  # The size of the double representation.
  set(NR_DOUBLE_SIZE ${DoubleSize})
endif()
message(STATUS "using double type: ${DoubleType}; size: ${DoubleSize}")

# defines used in qucs_typedefs.h
set(QUCS_INT32_TYPE ${nr_int32_t})
set(QUCS_INT16_TYPE ${nr_int16_t})
set(QUCS_DOUBLE_TYPE ${DoubleType})
set(QUCS_DOUBLE_SIZE ${DoubleSize})
