#ifndef _OURDEF_INCLUDED
#ifndef _SWEODEF_INCLUDED
#define _SWEODEF_INCLUDED

# define MY_TRUE 1
# define MY_FALSE 0

#ifdef __CYGWIN__
# undef __GNUC__
#endif

#if !defined(TLSOFF) && !defined( __APPLE__ ) && !defined(WIN32) && !defined(DOS32)
#if defined( __GNUC__ ) || defined( __CYGWIN__ )
#define TLS     __thread
#else
#define TLS     __declspec(thread)
#endif
#else
#define TLS
#endif

#ifdef _WIN32
# undef MSDOS
# define MSDOS MY_TRUE
#include <wtypes.h>
#include <objbase.h>
#include <wincon.h>
#include <winbase.h>
#include <io.h>
#include <windows.h>
# define sleep(x)	Sleep((x) * 1000)
#endif

#ifdef _MSC_VER
# define MS_VC
#endif

#ifdef WIN32
# define MSDOS MY_TRUE
#endif

#ifdef MSDOS
# undef MSDOS
# define MSDOS MY_TRUE
#endif

#ifdef __TURBOC__
# ifndef MSDOS
#   define MSDOS MY_TRUE
# endif
# define TURBO_C
#endif

#ifdef __SC__
# ifndef MSDOS
#   define MSDOS MY_TRUE
# endif
# define SYMANTEC_C
#endif

#ifdef __WATCOMC__
# ifndef MSDOS
#   define MSDOS MY_TRUE
# endif
# define WATCOMC
#endif

#ifdef MSDOS
#  define HPUNIX MY_FALSE
#  define INTEL_BYTE_ORDER 1
#  ifndef TURBO_C
#    define MS_C
#  endif
# define UNIX_FS MY_FALSE
#else
#  define MSDOS MY_FALSE
#  define HPUNIX MY_TRUE
#  ifndef _HPUX_SOURCE
#    define _HPUX_SOURCE
#  endif
#  define UNIX_FS MY_TRUE
#endif

#include <math.h>
#include <stdlib.h>
#ifndef FILE
# include <stdio.h>
#endif

#if HPUNIX
#  include <unistd.h>
#endif

#include <limits.h>
#if INT_MAX < 40000
# define INT_16
#else
# if LONG_MAX > INT_MAX
#   define LONG_64
# endif
#endif

#ifdef BYTE_ORDER
#ifdef LITTLE_ENDIAN
# if BYTE_ORDER == LITTLE_ENDIAN
#  define INTEL_BYTE_ORDER
# endif
#endif
#endif

#ifdef INT_16
  typedef long	int32;
  typedef unsigned long	uint32;
  typedef int	int16;
  typedef double  REAL8;
  typedef long    INT4;
  typedef unsigned long UINT4;

  typedef int     AS_BOOL;
  typedef unsigned int UINT2;
# define ABS4	labs
#else
  typedef int	int32;
  typedef long long	int64;
  typedef unsigned int	uint32;
  typedef short	int16;
  typedef double  REAL8;
  typedef int     INT4;
  typedef unsigned int UINT4;

  typedef int     AS_BOOL;
  typedef unsigned short UINT2;
  # define ABS4	abs
#endif

#if MSDOS
# ifdef TURBO_C
#   include <alloc.h>
# else
#   include <malloc.h>
# endif
# define SIGALRM SIGINT
#endif

#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#ifndef OK
#  define OK (0)
#  define ERR (-1)
#endif

#ifdef __GNUC__
#ifdef _WIN32
#define UCHAR SWE_UCHAR
#endif
#endif

typedef unsigned char UCHAR;
#define UCP	(UCHAR*)
#define SCP	(char*)

# define ODEGREE_STRING "°"

#ifndef HUGE
#  define HUGE 1.7E+308
#endif
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#define AS_MAXCH 256

#define RADTODEG (180.0 / M_PI)
#define DEGTORAD (M_PI / 180.0)

typedef int32    centisec;
#define CS	(centisec)
#define CSEC	centisec

#define DEG     360000
#define DEG7_30 (2700000)
#define DEG15   (15 * DEG)
#define DEG24   (24 * DEG)
#define DEG30   (30 * DEG)
#define DEG60   (60 * DEG)
#define DEG90   (90 * DEG)
#define DEG120  (120 * DEG)
#define DEG150  (150 * DEG)
#define DEG180  (180 * DEG)
#define DEG270  (270 * DEG)
#define DEG360  (360 * DEG)

#define CSTORAD	(DEGTORAD / 360000.0)
#define RADTOCS (RADTODEG * 360000.0)

#define CS2DEG	(1.0/360000.0)

#if UNIX_FS
#  define BFILE_R_ACCESS "r"
#  define BFILE_RW_ACCESS "r+"
#  define BFILE_W_CREATE "w"
#  define BFILE_A_ACCESS "a+"
#  define FILE_R_ACCESS "r"
#  define FILE_RW_ACCESS "r+"
#  define FILE_W_CREATE "w"
#  define FILE_A_ACCESS "a+"
#  define O_BINARY 0
#  define OPEN_MODE 0666
#  define DIR_GLUE "/"
#  define PATH_SEPARATOR ";:"
#else
#  define BFILE_R_ACCESS "rb"
#  define BFILE_RW_ACCESS "r+b"
#  define BFILE_W_CREATE "wb"
#  define BFILE_A_ACCESS "a+b"
#  define PATH_SEPARATOR ";"
#  define OPEN_MODE 0666
#  define FILE_R_ACCESS "rt"
#  define FILE_RW_ACCESS "r+t"
#  define FILE_W_CREATE "wt"
#  define FILE_A_ACCESS "a+t"

#  define DIR_GLUE "\\"
#endif

#include <string.h>
#include <ctype.h>

#endif
#endif
