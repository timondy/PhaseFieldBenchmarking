/* src/H5config.h.  Generated from H5config.h.in by configure.  */
/* src/H5config.h.in.  Generated from configure.in by autoheader.  */

/* Define if the memory buffers being written to disk should be cleared before
   writing. */
#define CLEAR_MEMORY 1

/* Define if your system can handle converting denormalized floating-point
   values. */
#define CONVERT_DENORMAL_FLOAT 1

/* Define if `dev_t' is a scalar */
#define DEV_T_IS_SCALAR 1

/* Define if gettimeofday() populates the tz pointer passed in */
#define GETTIMEOFDAY_GIVES_TZ 1

/* Define to 1 if you have the `alarm' function. */
#define HAVE_ALARM 1

/* Define if the __attribute__(()) extension is present */
#define HAVE_ATTRIBUTE 1

/* Define to 1 if you have the `BSDgettimeofday' function. */
/* #undef HAVE_BSDGETTIMEOFDAY */

/* Define to 1 if you have the declaration of `tzname', and to 0 if you don't.
   */
/* #undef HAVE_DECL_TZNAME */

/* Define to 1 if you have the `difftime' function. */
#define HAVE_DIFFTIME 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <dmalloc.h> header file. */
/* #undef HAVE_DMALLOC_H */

/* Define if library information should be embedded in the executables */
#define HAVE_EMBEDDED_LIBINFO 1

/* Define to 1 if you have the <features.h> header file. */
#define HAVE_FEATURES_H 1

/* Define if support for deflate filter is enabled */
/* #undef HAVE_FILTER_DEFLATE */

/* Define if support for Fletcher32 checksum is enabled */
#define HAVE_FILTER_FLETCHER32 1

/* Define if support for shuffle filter is enabled */
#define HAVE_FILTER_SHUFFLE 1

/* Define if support for szip filter is enabled */
/* #undef HAVE_FILTER_SZIP */

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Define to 1 if you have the `frexpf' function. */
#define HAVE_FREXPF 1

/* Define to 1 if you have the `frexpl' function. */
#define HAVE_FREXPL 1

/* Define to 1 if you have the `fseek64' function. */
/* #undef HAVE_FSEEK64 */

/* Define to 1 if you have the `fseeko' function. */
#define HAVE_FSEEKO 1

/* Define to 1 if you have the `ftello' function. */
#define HAVE_FTELLO 1

/* Define if the function stack tracing code is to be compiled in */
/* #undef HAVE_FUNCSTACK */

/* Define if the compiler understand the __FUNCTION__ keyword */
#define HAVE_FUNCTION 1

/* Define to 1 if you have the `GetConsoleScreenBufferInfo' function. */
/* #undef HAVE_GETCONSOLESCREENBUFFERINFO */

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to 1 if you have the `getpwuid' function. */
#define HAVE_GETPWUID 1

/* Define to 1 if you have the `getrusage' function. */
#define HAVE_GETRUSAGE 1

/* Define to 1 if you have the `gettextinfo' function. */
/* #undef HAVE_GETTEXTINFO */

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the `get_fpc_csr' function. */
/* #undef HAVE_GET_FPC_CSR */

/* Define if we have GPFS support */
/* #undef HAVE_GPFS */

/* Define to 1 if you have the <gpfs.h> header file. */
/* #undef HAVE_GPFS_H */

/* Define if library will contain instrumentation to detect correct
   optimization operation */
/* #undef HAVE_INSTRUMENTED_LIBRARY */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `ioctl' function. */
#define HAVE_IOCTL 1

/* Define to 1 if you have the <io.h> header file. */
/* #undef HAVE_IO_H */

/* Define to 1 if you have the `dmalloc' library (-ldmalloc). */
/* #undef HAVE_LIBDMALLOC */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `mpe' library (-lmpe). */
/* #undef HAVE_LIBMPE */

/* Define to 1 if you have the `mpi' library (-lmpi). */
/* #undef HAVE_LIBMPI */

/* Define to 1 if you have the `mpich' library (-lmpich). */
/* #undef HAVE_LIBMPICH */

/* Define to 1 if you have the `mpio' library (-lmpio). */
/* #undef HAVE_LIBMPIO */

/* Define to 1 if you have the `nsl' library (-lnsl). */
/* #undef HAVE_LIBNSL */

/* Define to 1 if you have the `pthread' library (-lpthread). */
/* #undef HAVE_LIBPTHREAD */

/* Define to 1 if you have the `sz' library (-lsz). */
/* #undef HAVE_LIBSZ */

/* Define to 1 if you have the `z' library (-lz). */
/* #undef HAVE_LIBZ */

/* Define to 1 if you have the `longjmp' function. */
#define HAVE_LONGJMP 1

/* Define to 1 if you have the `lseek64' function. */
#define HAVE_LSEEK64 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if we have MPE support */
/* #undef HAVE_MPE */

/* Define to 1 if you have the <mpe.h> header file. */
/* #undef HAVE_MPE_H */

/* Define if MPI_File_get_size works correctly */
#define HAVE_MPI_GET_SIZE 1

/* Define if we have parallel support */
#define HAVE_PARALLEL 1

/* Define to 1 if you have the <pthread.h> header file. */
/* #undef HAVE_PTHREAD_H */

/* Define to 1 if you have the <setjmp.h> header file. */
#define HAVE_SETJMP_H 1

/* Define to 1 if you have the `setsysinfo' function. */
/* #undef HAVE_SETSYSINFO */

/* Define to 1 if you have the `sigaction' function. */
#define HAVE_SIGACTION 1

/* Define to 1 if you have the `signal' function. */
#define HAVE_SIGNAL 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define if `struct stat' has the `st_blocks' field */
#define HAVE_STAT_ST_BLOCKS 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define if `struct text_info' is defined */
/* #undef HAVE_STRUCT_TEXT_INFO */

/* Define if `struct timezone' is defined */
#define HAVE_STRUCT_TIMEZONE 1

/* Define to 1 if `tm_zone' is member of `struct tm'. */
#define HAVE_STRUCT_TM_TM_ZONE 1

/* Define if `struct videoconfig' is defined */
/* #undef HAVE_STRUCT_VIDEOCONFIG */

/* Define to 1 if you have the `system' function. */
#define HAVE_SYSTEM 1

/* Define to 1 if you have the <sys/fpu.h> header file. */
/* #undef HAVE_SYS_FPU_H */

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#define HAVE_SYS_IOCTL_H 1

/* Define to 1 if you have the <sys/proc.h> header file. */
/* #undef HAVE_SYS_PROC_H */

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/sysinfo.h> header file. */
/* #undef HAVE_SYS_SYSINFO_H */

/* Define to 1 if you have the <sys/timeb.h> header file. */
#define HAVE_SYS_TIMEB_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <szlib.h> header file. */
/* #undef HAVE_SZLIB_H */

/* Define if we have thread safe support */
/* #undef HAVE_THREADSAFE */

/* Define if `timezone' is a global variable */
#define HAVE_TIMEZONE 1

/* Define if the ioctl TIOCGETD is defined */
#define HAVE_TIOCGETD 1

/* Define if the ioctl TIOGWINSZ is defined */
#define HAVE_TIOCGWINSZ 1

/* Define if `tm_gmtoff' is a member of `struct tm' */
#define HAVE_TM_GMTOFF 1

/* Define to 1 if your `struct tm' has `tm_zone'. Deprecated, use
   `HAVE_STRUCT_TM_TM_ZONE' instead. */
#define HAVE_TM_ZONE 1

/* Define to 1 if you don't have `tm_zone' but do have the external array
   `tzname'. */
/* #undef HAVE_TZNAME */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vsnprintf' function. */
#define HAVE_VSNPRINTF 1

/* Define to 1 if you have the `waitpid' function. */
#define HAVE_WAITPID 1

/* Define to 1 if you have the <winsock.h> header file. */
/* #undef HAVE_WINSOCK_H */

/* Define to 1 if you have the <zlib.h> header file. */
/* #undef HAVE_ZLIB_H */

/* Define to 1 if you have the `_getvideoconfig' function. */
/* #undef HAVE__GETVIDEOCONFIG */

/* Define to 1 if you have the `_scrsize' function. */
/* #undef HAVE__SCRSIZE */

/* Define if `__tm_gmtoff' is a member of `struct tm' */
/* #undef HAVE___TM_GMTOFF */

/* Define if your system can handle complicated MPI derived datatype
   correctly. */
/* #undef MPI_COMPLEX_DERIVED_DATATYPE_WORKS */

/* Define if your system's `MPI_File_set_size' function works for files over
   2GB. */
#define MPI_FILE_SET_SIZE_BIG 1

/* Define if shared writing must be disabled (CodeWarrior only) */
/* #undef NO_SHARED_WRITING */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "help@hdfgroup.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "HDF5"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "HDF5 1.6.10"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "hdf5"

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.6.10"

/* Width for printf() for type `long long' or `__int64', use `ll' */
#define PRINTF_LL_WIDTH "l"

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `int16_t', as computed by sizeof. */
#define SIZEOF_INT16_T 2

/* The size of `int32_t', as computed by sizeof. */
#define SIZEOF_INT32_T 4

/* The size of `int64_t', as computed by sizeof. */
#define SIZEOF_INT64_T 8

/* The size of `int8_t', as computed by sizeof. */
#define SIZEOF_INT8_T 1

/* The size of `int_fast16_t', as computed by sizeof. */
#define SIZEOF_INT_FAST16_T 8

/* The size of `int_fast32_t', as computed by sizeof. */
#define SIZEOF_INT_FAST32_T 8

/* The size of `int_fast64_t', as computed by sizeof. */
#define SIZEOF_INT_FAST64_T 8

/* The size of `int_fast8_t', as computed by sizeof. */
#define SIZEOF_INT_FAST8_T 1

/* The size of `int_least16_t', as computed by sizeof. */
#define SIZEOF_INT_LEAST16_T 2

/* The size of `int_least32_t', as computed by sizeof. */
#define SIZEOF_INT_LEAST32_T 4

/* The size of `int_least64_t', as computed by sizeof. */
#define SIZEOF_INT_LEAST64_T 8

/* The size of `int_least8_t', as computed by sizeof. */
#define SIZEOF_INT_LEAST8_T 1

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long double', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE 16

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `off_t', as computed by sizeof. */
#define SIZEOF_OFF_T 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* The size of `ssize_t', as computed by sizeof. */
#define SIZEOF_SSIZE_T 8

/* The size of `uint16_t', as computed by sizeof. */
#define SIZEOF_UINT16_T 2

/* The size of `uint32_t', as computed by sizeof. */
#define SIZEOF_UINT32_T 4

/* The size of `uint64_t', as computed by sizeof. */
#define SIZEOF_UINT64_T 8

/* The size of `uint8_t', as computed by sizeof. */
#define SIZEOF_UINT8_T 1

/* The size of `uint_fast16_t', as computed by sizeof. */
#define SIZEOF_UINT_FAST16_T 8

/* The size of `uint_fast32_t', as computed by sizeof. */
#define SIZEOF_UINT_FAST32_T 8

/* The size of `uint_fast64_t', as computed by sizeof. */
#define SIZEOF_UINT_FAST64_T 8

/* The size of `uint_fast8_t', as computed by sizeof. */
#define SIZEOF_UINT_FAST8_T 1

/* The size of `uint_least16_t', as computed by sizeof. */
#define SIZEOF_UINT_LEAST16_T 2

/* The size of `uint_least32_t', as computed by sizeof. */
#define SIZEOF_UINT_LEAST32_T 4

/* The size of `uint_least64_t', as computed by sizeof. */
#define SIZEOF_UINT_LEAST64_T 8

/* The size of `uint_least8_t', as computed by sizeof. */
#define SIZEOF_UINT_LEAST8_T 1

/* The size of `__int64', as computed by sizeof. */
#define SIZEOF___INT64 0

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if strict file format checks are enabled */
/* #undef STRICT_FORMAT_CHECKS */

/* Define if your system supports pthread_attr_setscope(&attribute,
   PTHREAD_SCOPE_SYSTEM) call. */
#define SYSTEM_SCOPE_THREADS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* Define if the HDF5 v1.4 compatibility functions are to be compiled in */
/* #undef WANT_H5_V1_4_COMPAT */

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#define inline __inline__
#endif

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `unsigned long' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `long' if <sys/types.h> does not define. */
/* #undef ssize_t */

#if defined(__cplusplus) && defined(inline)
#undef inline
#endif
