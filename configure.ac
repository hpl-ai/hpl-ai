AC_PREREQ([2.69])

AC_INIT(hpl-ai, 0.1, hpl-ai@icl.utk.edu)
AC_CONFIG_SRCDIR([hpl-ai.h])
AC_CONFIG_HEADERS([hpl-ai-config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_INSTALL

AM_INIT_AUTOMAKE([subdir-objects])

AM_PROG_CC_STDC
AM_PROG_CC_C_O

AC_CHECK_HEADERS([sys/time.h math.h])

AC_SEARCH_LIBS(sqrt, m, , AC_MSG_ERROR([sqrt() function not found in standard C math library]))

AC_CHECK_FUNCS([fabs])

AC_CHECK_HEADER([cuda.h], [], AC_MSG_FAILURE([Couldn't find cuda.h]))
AC_CHECK_HEADER([cuda_runtime.h], [], AC_MSG_FAILURE([Couldn't find cuda_runtime.h]))
AC_CHECK_HEADER([cublas_v2.h], [], AC_MSG_FAILURE([Couldn't find cublas_v2.h]))
AC_CHECK_LIB([cublas], [cublasAlloc], [], AC_MSG_FAILURE([Couldn't find libcublas]))
AC_CHECK_LIB([cudart], [cudaSetDevice], [], AC_MSG_FAILURE([Couldn't find libcudart]))

AC_CONFIG_FILES([Makefile])


AC_OUTPUT
