AU_ALIAS([ACX_SIRIUS], [AX_SIRIUS])
AC_DEFUN([ACX_SIRIUS],
         [
         AC_ARG_WITH(sirius,
                     [AC_HELP_STRING([--with-sirius=yes/no/prefix], [enable sirius planewave dft support (default yes)])],
                     [sirius_prefix="$withval"],
                     [sirius_prefix="/usr/local"])
         sirius_found=no
         sirius_headers_found=no
         AS_IF([test "$with_sirius" != "no" -a "${MPICC}yes" != "yes"],
               [
               AC_CHECK_FILE(${sirius_prefix}/include/sirius/sirius.h, [sirius_headers_found=yes], [sirius_headers_found=no])
               AC_CHECK_FILE(${sirius_prefix}/lib64/libsirius.a, [sirius_libs_found=yes], [sirius_libs_found=no])
#               AC_CHECK_FILE(${sirius_prefix}/lib64/libsirius_cu.a,
#                              [SIRIUS_CUDA_LDFLAGS="-L${sirius_prefix}/lib"
#                              SIRIUS_CUDA_LIBS="-lsirius"
#                              sirius_cuda_libs_found=yes],
#                              [sirius_cuda_libs_found=no])
               AS_IF([test "$sirius_libs_found" == "yes" -a "$sirius_headers_found" == "yes"],
                     [SIRIUS_CFLAGS="-I${sirius_prefix}/include"
                      SIRIUS_LDFLAGS="-L${sirius_prefix}/lib64"
                      SIRIUS_LIBS="-lsirius"
                     ])
               ])
         AC_SUBST(SIRIUS_CFLAGS)
         AC_SUBST(SIRIUS_LIBS)
         AC_SUBST(SIRIUS_LDFLAGS)
         AS_IF([test "$sirius_libs_found" == "yes" -a "$sirius_headers_found" == "yes"],
               [$1],
               [$2])
])
