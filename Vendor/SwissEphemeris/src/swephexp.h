#ifdef __cplusplus
extern "C" {
#endif

#ifndef _SWEPHEXP_INCLUDED
#define _SWEPHEXP_INCLUDED

#include "sweodef.h"

#define SE_AUNIT_TO_KM        (149597870.700)
#define SE_AUNIT_TO_LIGHTYEAR (1.0/63241.07708427)
#define SE_AUNIT_TO_PARSEC    (1.0/206264.8062471)

# define SE_JUL_CAL	0
# define SE_GREG_CAL	1

#define SE_ECL_NUT      -1

#define SE_SUN          0
#define SE_MOON         1
#define SE_MERCURY      2
#define SE_VENUS        3
#define SE_MARS         4
#define SE_JUPITER      5
#define SE_SATURN       6
#define SE_URANUS       7
#define SE_NEPTUNE      8
#define SE_PLUTO        9
#define SE_MEAN_NODE    10
#define SE_TRUE_NODE    11
#define SE_MEAN_APOG    12
#define SE_OSCU_APOG    13
#define SE_EARTH        14
#define SE_CHIRON       15
#define SE_PHOLUS       16
#define SE_CERES        17
#define SE_PALLAS       18
#define SE_JUNO         19
#define SE_VESTA        20
#define SE_INTP_APOG    21
#define SE_INTP_PERG    22

#define SE_NPLANETS     23

#define SE_PLMOON_OFFSET   9000
#define SE_AST_OFFSET   10000
#define SE_VARUNA   (SE_AST_OFFSET + 20000)

#define SE_FICT_OFFSET  	40
#define SE_FICT_OFFSET_1  	39
#define SE_FICT_MAX  	       999
#define SE_NFICT_ELEM           15

#define SE_COMET_OFFSET 1000

#define SE_NALL_NAT_POINTS      (SE_NPLANETS + SE_NFICT_ELEM)

#define SE_CUPIDO       	40
#define SE_HADES        	41
#define SE_ZEUS         	42
#define SE_KRONOS       	43
#define SE_APOLLON      	44
#define SE_ADMETOS      	45
#define SE_VULKANUS     	46
#define SE_POSEIDON     	47

#define SE_ISIS         	48
#define SE_NIBIRU       	49
#define SE_HARRINGTON           50
#define SE_NEPTUNE_LEVERRIER    51
#define SE_NEPTUNE_ADAMS        52
#define SE_PLUTO_LOWELL         53
#define SE_PLUTO_PICKERING      54
#define SE_VULCAN      		55
#define SE_WHITE_MOON  		56
#define SE_PROSERPINA  		57
#define SE_WALDEMATH  		58

#define SE_FIXSTAR      -10

#define SE_ASC			0
#define SE_MC			1
#define SE_ARMC			2
#define SE_VERTEX		3
#define SE_EQUASC  		4
#define SE_COASC1		5
#define SE_COASC2		6
#define SE_POLASC		7
#define SE_NASCMC		8

#define SEFLG_JPLEPH    1
#define SEFLG_SWIEPH    2
#define SEFLG_MOSEPH    4

#define SEFLG_HELCTR	8
#define SEFLG_TRUEPOS	16
#define SEFLG_J2000	32
#define SEFLG_NONUT	64
#define SEFLG_SPEED3	128

#define SEFLG_SPEED	256
#define SEFLG_NOGDEFL	512
#define SEFLG_NOABERR	1024
#define SEFLG_ASTROMETRIC (SEFLG_NOABERR|SEFLG_NOGDEFL)

#define SEFLG_EQUATORIAL (2*1024)
#define SEFLG_XYZ	(4*1024)
#define SEFLG_RADIANS	(8*1024)
#define SEFLG_BARYCTR	(16*1024)
#define SEFLG_TOPOCTR	(32*1024)
#define SEFLG_ORBEL_AA SEFLG_TOPOCTR

#define SEFLG_TROPICAL	(0)
#define SEFLG_SIDEREAL	(64*1024)
#define SEFLG_ICRS	(128*1024)
#define SEFLG_DPSIDEPS_1980	(256*1024)

#define SEFLG_JPLHOR	SEFLG_DPSIDEPS_1980
#define SEFLG_JPLHOR_APPROX	(512*1024)
#define SEFLG_CENTER_BODY	(1024*1024)

#define SEFLG_TEST_PLMOON	(2*1024*1024 | SEFLG_J2000 | SEFLG_ICRS | SEFLG_HELCTR | SEFLG_TRUEPOS)

#define SE_SIDBITS		256

#define SE_SIDBIT_ECL_T0        256

#define SE_SIDBIT_SSY_PLANE     512

#define SE_SIDBIT_USER_UT       1024

#define SE_SIDBIT_ECL_DATE      2048

#define SE_SIDBIT_NO_PREC_OFFSET       4096

#define SE_SIDBIT_PREC_ORIG     8192

#define SE_SIDM_FAGAN_BRADLEY    0
#define SE_SIDM_LAHIRI           1
#define SE_SIDM_DELUCE           2
#define SE_SIDM_RAMAN            3
#define SE_SIDM_USHASHASHI       4
#define SE_SIDM_KRISHNAMURTI     5
#define SE_SIDM_DJWHAL_KHUL      6
#define SE_SIDM_YUKTESHWAR       7
#define SE_SIDM_JN_BHASIN        8
#define SE_SIDM_BABYL_KUGLER1    9
#define SE_SIDM_BABYL_KUGLER2   10
#define SE_SIDM_BABYL_KUGLER3   11
#define SE_SIDM_BABYL_HUBER    	12
#define SE_SIDM_BABYL_ETPSC    	13
#define SE_SIDM_ALDEBARAN_15TAU 14
#define SE_SIDM_HIPPARCHOS      15
#define SE_SIDM_SASSANIAN       16
#define SE_SIDM_GALCENT_0SAG    17
#define SE_SIDM_J2000           18
#define SE_SIDM_J1900           19
#define SE_SIDM_B1950           20
#define SE_SIDM_SURYASIDDHANTA  21
#define SE_SIDM_SURYASIDDHANTA_MSUN  22
#define SE_SIDM_ARYABHATA       23
#define SE_SIDM_ARYABHATA_MSUN  24
#define SE_SIDM_SS_REVATI       25
#define SE_SIDM_SS_CITRA        26
#define SE_SIDM_TRUE_CITRA      27
#define SE_SIDM_TRUE_REVATI     28
#define SE_SIDM_TRUE_PUSHYA     29
#define SE_SIDM_GALCENT_RGILBRAND 30
#define SE_SIDM_GALEQU_IAU1958  31
#define SE_SIDM_GALEQU_TRUE     32
#define SE_SIDM_GALEQU_MULA     33
#define SE_SIDM_GALALIGN_MARDYKS 34
#define SE_SIDM_TRUE_MULA       35
#define SE_SIDM_GALCENT_MULA_WILHELM       36
#define SE_SIDM_ARYABHATA_522   37
#define SE_SIDM_BABYL_BRITTON   38
#define SE_SIDM_TRUE_SHEORAN  	39
#define SE_SIDM_GALCENT_COCHRANE   	40
#define SE_SIDM_GALEQU_FIORENZA 41
#define SE_SIDM_VALENS_MOON     42
#define SE_SIDM_LAHIRI_1940     43
#define SE_SIDM_LAHIRI_VP285    44
#define SE_SIDM_KRISHNAMURTI_VP291    45
#define SE_SIDM_LAHIRI_ICRC     46

#define SE_SIDM_USER            255

#define SE_NSIDM_PREDEF	        47

#define SE_NODBIT_MEAN		1
#define SE_NODBIT_OSCU		2
#define SE_NODBIT_OSCU_BAR	4
#define SE_NODBIT_FOPOINT	256

#define SEFLG_DEFAULTEPH SEFLG_SWIEPH

#define SE_MAX_STNAME		256

#define SE_ECL_CENTRAL		1
#define SE_ECL_NONCENTRAL	2
#define SE_ECL_TOTAL		4
#define SE_ECL_ANNULAR		8
#define SE_ECL_PARTIAL		16
#define SE_ECL_ANNULAR_TOTAL	32
#define SE_ECL_HYBRID   	32
#define SE_ECL_PENUMBRAL	64
#define SE_ECL_ALLTYPES_SOLAR   (SE_ECL_CENTRAL|SE_ECL_NONCENTRAL|SE_ECL_TOTAL|SE_ECL_ANNULAR|SE_ECL_PARTIAL|SE_ECL_ANNULAR_TOTAL)
#define SE_ECL_ALLTYPES_LUNAR   (SE_ECL_TOTAL|SE_ECL_PARTIAL|SE_ECL_PENUMBRAL)
#define SE_ECL_VISIBLE			128
#define SE_ECL_MAX_VISIBLE		256
#define SE_ECL_1ST_VISIBLE		512
#define SE_ECL_PARTBEG_VISIBLE		512
#define SE_ECL_2ND_VISIBLE		1024
#define SE_ECL_TOTBEG_VISIBLE		1024
#define SE_ECL_3RD_VISIBLE		2048
#define SE_ECL_TOTEND_VISIBLE		2048
#define SE_ECL_4TH_VISIBLE		4096
#define SE_ECL_PARTEND_VISIBLE		4096
#define SE_ECL_PENUMBBEG_VISIBLE	8192
#define SE_ECL_PENUMBEND_VISIBLE	16384
#define SE_ECL_OCC_BEG_DAYLIGHT		8192
#define SE_ECL_OCC_END_DAYLIGHT		16384
#define SE_ECL_ONE_TRY          (32*1024)

#define SE_CALC_RISE		1
#define SE_CALC_SET		2
#define SE_CALC_MTRANSIT	4
#define SE_CALC_ITRANSIT	8
#define SE_BIT_DISC_CENTER      256

#define SE_BIT_DISC_BOTTOM      8192

#define SE_BIT_GEOCTR_NO_ECL_LAT 128

#define SE_BIT_NO_REFRACTION    512

#define SE_BIT_CIVIL_TWILIGHT    1024
#define SE_BIT_NAUTIC_TWILIGHT   2048
#define SE_BIT_ASTRO_TWILIGHT    4096
#define SE_BIT_FIXED_DISC_SIZE  16384

#define SE_BIT_FORCE_SLOW_METHOD 32768

#define SE_BIT_HINDU_RISING  (SE_BIT_DISC_CENTER|SE_BIT_NO_REFRACTION|SE_BIT_GEOCTR_NO_ECL_LAT)

#define SE_ECL2HOR		0
#define SE_EQU2HOR		1
#define SE_HOR2ECL		0
#define SE_HOR2EQU		1

#define SE_TRUE_TO_APP	0
#define SE_APP_TO_TRUE	1

#define SE_DE_NUMBER    431
#define SE_FNAME_DE200  "de200.eph"
#define SE_FNAME_DE403  "de403.eph"
#define SE_FNAME_DE404  "de404.eph"
#define SE_FNAME_DE405  "de405.eph"
#define SE_FNAME_DE406  "de406.eph"
#define SE_FNAME_DE431  "de431.eph"
#define SE_FNAME_DFT    SE_FNAME_DE431
#define SE_FNAME_DFT2   SE_FNAME_DE406
#define SE_STARFILE_OLD "fixstars.cat"
#define SE_STARFILE     "sefstars.txt"
#define SE_ASTNAMFILE   "seasnam.txt"
#define SE_FICTFILE     "seorbel.txt"

#ifndef SE_EPHE_PATH
#if MSDOS
#  define SE_EPHE_PATH    "\\sweph\\ephe\\"
#else
#  define SE_EPHE_PATH    ".:/users/ephe2/:/users/ephe/"

#endif
#endif

# define SE_SPLIT_DEG_ROUND_SEC    1
# define SE_SPLIT_DEG_ROUND_MIN    2
# define SE_SPLIT_DEG_ROUND_DEG    4
# define SE_SPLIT_DEG_ZODIACAL     8
# define SE_SPLIT_DEG_NAKSHATRA 1024
# define SE_SPLIT_DEG_KEEP_SIGN   16

# define SE_SPLIT_DEG_KEEP_DEG    32

#define SE_HELIACAL_RISING		1
#define SE_HELIACAL_SETTING		2
#define SE_MORNING_FIRST		SE_HELIACAL_RISING
#define SE_EVENING_LAST			SE_HELIACAL_SETTING
#define SE_EVENING_FIRST		3
#define SE_MORNING_LAST			4
#define SE_ACRONYCHAL_RISING		5
#define SE_ACRONYCHAL_SETTING		6
#define SE_COSMICAL_SETTING		SE_ACRONYCHAL_SETTING

#define SE_HELFLAG_LONG_SEARCH 		128
#define SE_HELFLAG_HIGH_PRECISION 	256
#define SE_HELFLAG_OPTICAL_PARAMS	512
#define SE_HELFLAG_NO_DETAILS		1024
#define SE_HELFLAG_SEARCH_1_PERIOD	(1 << 11)
#define SE_HELFLAG_VISLIM_DARK		(1 << 12)
#define SE_HELFLAG_VISLIM_NOMOON	(1 << 13)

#define SE_HELFLAG_VISLIM_PHOTOPIC	(1 << 14)
#define SE_HELFLAG_VISLIM_SCOTOPIC	(1 << 15)
#define SE_HELFLAG_AV	 		(1 << 16)
#define SE_HELFLAG_AVKIND_VR 		(1 << 16)
#define SE_HELFLAG_AVKIND_PTO 		(1 << 17)
#define SE_HELFLAG_AVKIND_MIN7 		(1 << 18)
#define SE_HELFLAG_AVKIND_MIN9 		(1 << 19)
#define SE_HELFLAG_AVKIND (SE_HELFLAG_AVKIND_VR|SE_HELFLAG_AVKIND_PTO|SE_HELFLAG_AVKIND_MIN7|SE_HELFLAG_AVKIND_MIN9)
#define TJD_INVALID		 	99999999.0
#define SIMULATE_VICTORVB               1

#if 0
#define SE_HELIACAL_LONG_SEARCH 	128
#define SE_HELIACAL_HIGH_PRECISION 	256
#define SE_HELIACAL_OPTICAL_PARAMS	512
#define SE_HELIACAL_NO_DETAILS		1024
#define SE_HELIACAL_SEARCH_1_PERIOD	(1 << 11)
#define SE_HELIACAL_VISLIM_DARK		(1 << 12)
#define SE_HELIACAL_VISLIM_NOMOON	(1 << 13)
#define SE_HELIACAL_VISLIM_PHOTOPIC	(1 << 14)
#define SE_HELIACAL_AVKIND_VR 		(1 << 15)
#define SE_HELIACAL_AVKIND_PTO 		(1 << 16)
#define SE_HELIACAL_AVKIND_MIN7 		(1 << 17)
#define SE_HELIACAL_AVKIND_MIN9 		(1 << 18)
#define SE_HELIACAL_AVKIND (SE_HELFLAG_AVKIND_VR|SE_HELFLAG_AVKIND_PTO|SE_HELFLAG_AVKIND_MIN7|SE_HELFLAG_AVKIND_MIN9)
#endif

#define SE_PHOTOPIC_FLAG		0
#define SE_SCOTOPIC_FLAG		1
#define SE_MIXEDOPIC_FLAG		2

#define SE_TIDAL_DE200          (-23.8946)
#define SE_TIDAL_DE403          (-25.580)
#define SE_TIDAL_DE404          (-25.580)
#define SE_TIDAL_DE405          (-25.826)
#define SE_TIDAL_DE406          (-25.826)
#define SE_TIDAL_DE421          (-25.85)
#define SE_TIDAL_DE422          (-25.85)
#define SE_TIDAL_DE430          (-25.82)
#define SE_TIDAL_DE431          (-25.80)
#define SE_TIDAL_DE441          (-25.936)
#define SE_TIDAL_26             (-26.0)
#define SE_TIDAL_STEPHENSON_2016             (-25.85)
#define SE_TIDAL_DEFAULT        SE_TIDAL_DE431
#define SE_TIDAL_AUTOMATIC             999999
#define SE_TIDAL_MOSEPH                SE_TIDAL_DE404
#define SE_TIDAL_SWIEPH                SE_TIDAL_DEFAULT
#define SE_TIDAL_JPLEPH                SE_TIDAL_DEFAULT

#define SE_DELTAT_AUTOMATIC             (-1E-10)

#define SE_MODEL_DELTAT         0
#define SE_MODEL_PREC_LONGTERM  1
#define SE_MODEL_PREC_SHORTTERM 2
#define SE_MODEL_NUT            3
#define SE_MODEL_BIAS           4
#define SE_MODEL_JPLHOR_MODE    5
#define SE_MODEL_JPLHORA_MODE   6
#define SE_MODEL_SIDT           7
#define NSE_MODELS              8

#define SEMOD_NPREC		11
#define SEMOD_PREC_IAU_1976      1
#define SEMOD_PREC_LASKAR_1986   2
#define SEMOD_PREC_WILL_EPS_LASK 3
#define SEMOD_PREC_WILLIAMS_1994 4
#define SEMOD_PREC_SIMON_1994    5
#define SEMOD_PREC_IAU_2000      6
#define SEMOD_PREC_BRETAGNON_2003      7
#define SEMOD_PREC_IAU_2006      8
#define SEMOD_PREC_VONDRAK_2011  9
#define SEMOD_PREC_OWEN_1990     10
#define SEMOD_PREC_NEWCOMB       11
#define SEMOD_PREC_DEFAULT       SEMOD_PREC_VONDRAK_2011

#define SEMOD_PREC_DEFAULT_SHORT SEMOD_PREC_VONDRAK_2011

#define SEMOD_NNUT		5
#define SEMOD_NUT_IAU_1980          1
#define SEMOD_NUT_IAU_CORR_1987     2

#define SEMOD_NUT_IAU_2000A         3
#define SEMOD_NUT_IAU_2000B         4
#define SEMOD_NUT_WOOLARD           5
#define SEMOD_NUT_DEFAULT           SEMOD_NUT_IAU_2000B

#define SEMOD_NSIDT		4
#define SEMOD_SIDT_IAU_1976         1
#define SEMOD_SIDT_IAU_2006         2
#define SEMOD_SIDT_IERS_CONV_2010   3
#define SEMOD_SIDT_LONGTERM         4
#define SEMOD_SIDT_DEFAULT          SEMOD_SIDT_LONGTERM

#define SEMOD_NBIAS		3
#define SEMOD_BIAS_NONE             1
#define SEMOD_BIAS_IAU2000          2
#define SEMOD_BIAS_IAU2006          3
#define SEMOD_BIAS_DEFAULT          SEMOD_BIAS_IAU2006

#define SEMOD_NJPLHOR		2
#define SEMOD_JPLHOR_LONG_AGREEMENT  1

#define SEMOD_JPLHOR_DEFAULT        SEMOD_JPLHOR_LONG_AGREEMENT

#define SEMOD_NJPLHORA		3
#define SEMOD_JPLHORA_1     1
#define SEMOD_JPLHORA_2     2
#define SEMOD_JPLHORA_3     3
#define SEMOD_JPLHORA_DEFAULT     SEMOD_JPLHORA_3

#define SEMOD_NDELTAT		5
#define SEMOD_DELTAT_STEPHENSON_MORRISON_1984   1
#define SEMOD_DELTAT_STEPHENSON_1997   2
#define SEMOD_DELTAT_STEPHENSON_MORRISON_2004   3
#define SEMOD_DELTAT_ESPENAK_MEEUS_2006   4
#define SEMOD_DELTAT_STEPHENSON_ETC_2016   5

#define SEMOD_DELTAT_DEFAULT   SEMOD_DELTAT_STEPHENSON_ETC_2016

#if defined(MAKE_DLL) || defined(USE_DLL) || defined(_WINDOWS)
#  include <windows.h>
extern HANDLE dllhandle;

#endif

#ifdef USE_DLL
#  include "swedll.h"
#endif

#if defined(DOS32) || !MSDOS || defined(WIN32)

#  define MALLOC malloc
#  define CALLOC calloc
#  define FREE free
#else
#  ifdef __BORLANDC__
#    include <alloc.h>
#    define MALLOC farmalloc
#    define CALLOC farcalloc
#    define FREE farfree
#  else
#    define MALLOC _fmalloc
#    define CALLOC _fcalloc
#    define FREE _ffree
#  endif
#endif

#ifdef MAKE_DLL
  #if defined (PASCAL) || defined(__stdcall)
   #if defined UNDECO_DLL
    #define CALL_CONV __cdecl
   #else
    #define CALL_CONV __stdcall
   #endif
  #else
    #define CALL_CONV
  #endif

  #define EXP32  __declspec( dllexport )
#else
  #define CALL_CONV
  #define EXP32
#endif

#ifndef _SWEDLL_H

#define ext_def(x)	extern EXP32 x CALL_CONV

ext_def(int32) swe_heliacal_ut(double tjdstart_ut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 iflag, double *dret, char *serr);
ext_def(int32) swe_heliacal_pheno_ut(double tjd_ut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 helflag, double *darr, char *serr);
ext_def(int32) swe_vis_limit_mag(double tjdut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 helflag, double *dret, char *serr);

ext_def(int32) swe_heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr);
ext_def(int32) swe_topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr);

ext_def(void) swe_set_astro_models(char *samod, int32 iflag);
ext_def(void) swe_get_astro_models(char *samod, char *sdet, int32 iflag);

ext_def(char *) swe_version(char *);
ext_def(char *) swe_get_library_path(char *);

ext_def( int32 ) swe_calc(
        double tjd, int ipl, int32 iflag,
        double *xx,
        char *serr);

ext_def(int32) swe_calc_ut(double tjd_ut, int32 ipl, int32 iflag,
	double *xx, char *serr);

ext_def(int32) swe_calc_pctr(double tjd, int32 ipl, int32 iplctr, int32 iflag, double *xxret, char *serr);

ext_def(double) swe_solcross(double x2cross, double jd_et, int32 flag, char *serr);
ext_def(double) swe_solcross_ut(double x2cross, double jd_ut, int32 flag, char *serr);
ext_def(double) swe_mooncross(double x2cross, double jd_et, int32 flag, char *serr);
ext_def(double) swe_mooncross_ut(double x2cross, double jd_ut, int32 flag, char *serr);
ext_def(double) swe_mooncross_node(double jd_et, int32 flag, double *xlon, double *xlat, char *serr);
ext_def(double) swe_mooncross_node_ut(double jd_ut, int32 flag, double *xlon, double *xlat, char *serr);
ext_def(int32) swe_helio_cross(int32 ipl, double x2cross, double jd_et, int32 iflag, int32 dir, double *jd_cross, char *serr);
ext_def(int32) swe_helio_cross_ut(int32 ipl, double x2cross, double jd_ut, int32 iflag, int32 dir, double *jd_cross, char *serr);

ext_def( int32 ) swe_fixstar(
        char *star, double tjd, int32 iflag,
        double *xx,
        char *serr);

ext_def(int32) swe_fixstar_ut(char *star, double tjd_ut, int32 iflag,
	double *xx, char *serr);

ext_def(int32) swe_fixstar_mag(char *star, double *mag, char *serr);

ext_def( int32 ) swe_fixstar2(
        char *star, double tjd, int32 iflag,
        double *xx,
        char *serr);

ext_def(int32) swe_fixstar2_ut(char *star, double tjd_ut, int32 iflag,
	double *xx, char *serr);

ext_def(int32) swe_fixstar2_mag(char *star, double *mag, char *serr);

ext_def( void ) swe_close(void);

ext_def( void ) swe_set_ephe_path(const char *path);

ext_def( void ) swe_set_jpl_file(const char *fname);

ext_def( char *) swe_get_planet_name(int ipl, char *spname);

ext_def (void) swe_set_topo(double geolon, double geolat, double geoalt);

ext_def(void) swe_set_sid_mode(int32 sid_mode, double t0, double ayan_t0);

ext_def(int32) swe_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr);
ext_def(int32) swe_get_ayanamsa_ex_ut(double tjd_ut, int32 iflag, double *daya, char *serr);
ext_def(double) swe_get_ayanamsa(double tjd_et);
ext_def(double) swe_get_ayanamsa_ut(double tjd_ut);

ext_def(const char *) swe_get_ayanamsa_name(int32 isidmode);
ext_def(const char *) swe_get_current_file_data(int ifno, double *tfstart, double *tfend, int *denum);

ext_def( int ) swe_date_conversion(
        int y , int m , int d ,
        double utime,
        char c,
        double *tjd);

ext_def( double ) swe_julday(
        int year, int month, int day, double hour,
        int gregflag);

ext_def( void ) swe_revjul (
        double jd,
        int gregflag,
        int *jyear, int *jmon, int *jday, double *jut);

ext_def(int32) swe_utc_to_jd(
        int32 iyear, int32 imonth, int32 iday,
	int32 ihour, int32 imin, double dsec,
	int32 gregflag, double *dret, char *serr);

ext_def(void) swe_jdet_to_utc(
        double tjd_et, int32 gregflag,
	int32 *iyear, int32 *imonth, int32 *iday,
	int32 *ihour, int32 *imin, double *dsec);

ext_def(void) swe_jdut1_to_utc(
        double tjd_ut, int32 gregflag,
	int32 *iyear, int32 *imonth, int32 *iday,
	int32 *ihour, int32 *imin, double *dsec);

ext_def(void) swe_utc_time_zone(
        int32 iyear, int32 imonth, int32 iday,
	int32 ihour, int32 imin, double dsec,
	double d_timezone,
	int32 *iyear_out, int32 *imonth_out, int32 *iday_out,
	int32 *ihour_out, int32 *imin_out, double *dsec_out);

ext_def( int ) swe_houses(
        double tjd_ut, double geolat, double geolon, int hsys,
	double *cusps, double *ascmc);

ext_def( int ) swe_houses_ex(
        double tjd_ut, int32 iflag, double geolat, double geolon, int hsys,
	double *cusps, double *ascmc);

ext_def( int ) swe_houses_ex2(
        double tjd_ut, int32 iflag, double geolat, double geolon, int hsys,
	double *cusps, double *ascmc, double *cusp_speed, double *ascmc_speed, char *serr);

ext_def( int ) swe_houses_armc(
        double armc, double geolat, double eps, int hsys,
	double *cusps, double *ascmc);

ext_def( int ) swe_houses_armc_ex2(
        double armc, double geolat, double eps, int hsys,
	double *cusps, double *ascmc, double *cusp_speed, double *ascmc_speed, char *serr);

ext_def(double) swe_house_pos(
	double armc, double geolat, double eps, int hsys, double *xpin, char *serr);

ext_def(const char *) swe_house_name(int hsys);

ext_def(int32) swe_gauquelin_sector(double t_ut, int32 ipl, char *starname, int32 iflag, int32 imeth, double *geopos, double atpress, double attemp, double *dgsect, char *serr);

ext_def (int32) swe_sol_eclipse_where(double tjd, int32 ifl, double *geopos, double *attr, char *serr);

ext_def (int32) swe_lun_occult_where(double tjd, int32 ipl, char *starname, int32 ifl, double *geopos, double *attr, char *serr);

ext_def (int32) swe_sol_eclipse_how(double tjd, int32 ifl, double *geopos, double *attr, char *serr);

ext_def (int32) swe_sol_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr);

ext_def (int32) swe_lun_occult_when_loc(double tjd_start, int32 ipl, char *starname, int32 ifl,
     double *geopos, double *tret, double *attr, int32 backward, char *serr);

ext_def (int32) swe_sol_eclipse_when_glob(double tjd_start, int32 ifl, int32 ifltype,
     double *tret, int32 backward, char *serr);

ext_def (int32) swe_lun_occult_when_glob(double tjd_start, int32 ipl, char *starname, int32 ifl, int32 ifltype,
     double *tret, int32 backward, char *serr);

ext_def (int32) swe_lun_eclipse_how(
          double tjd_ut,
          int32 ifl,
          double *geopos,
          double *attr,
          char *serr);

ext_def (int32) swe_lun_eclipse_when(double tjd_start, int32 ifl, int32 ifltype,
     double *tret, int32 backward, char *serr);

ext_def (int32) swe_lun_eclipse_when_loc(double tjd_start, int32 ifl,
     double *geopos, double *tret, double *attr, int32 backward, char *serr);

ext_def (int32) swe_pheno(double tjd, int32 ipl, int32 iflag, double *attr, char *serr);

ext_def(int32) swe_pheno_ut(double tjd_ut, int32 ipl, int32 iflag, double *attr, char *serr);

ext_def (double) swe_refrac(double inalt, double atpress, double attemp, int32 calc_flag);

ext_def (double) swe_refrac_extended(double inalt, double geoalt, double atpress, double attemp, double lapse_rate, int32 calc_flag, double *dret);

ext_def (void) swe_set_lapse_rate(double lapse_rate);

ext_def (void) swe_azalt(
      double tjd_ut,
      int32 calc_flag,
      double *geopos,
      double atpress,
      double attemp,
      double *xin,
      double *xaz);

ext_def (void) swe_azalt_rev(
      double tjd_ut,
      int32 calc_flag,
      double *geopos,
      double *xin,
      double *xout);

ext_def (int32) swe_rise_trans_true_hor(
               double tjd_ut, int32 ipl, char *starname,
	       int32 epheflag, int32 rsmi,
               double *geopos,
	       double atpress, double attemp,
	       double horhgt,
               double *tret,
               char *serr);

ext_def (int32) swe_rise_trans(
               double tjd_ut, int32 ipl, char *starname,
	       int32 epheflag, int32 rsmi,
               double *geopos,
	       double atpress, double attemp,
               double *tret,
               char *serr);

ext_def (int32) swe_nod_aps(double tjd_et, int32 ipl, int32 iflag,
                      int32  method,
                      double *xnasc, double *xndsc,
                      double *xperi, double *xaphe,
                      char *serr);

ext_def (int32) swe_nod_aps_ut(double tjd_ut, int32 ipl, int32 iflag,
                      int32  method,
                      double *xnasc, double *xndsc,
                      double *xperi, double *xaphe,
                      char *serr);
ext_def (int32) swe_get_orbital_elements(
  double tjd_et, int32 ipl, int32 iflag, double *dret, char *serr);

ext_def (int32) swe_orbit_max_min_true_distance(double tjd_et, int32 ipl, int32 iflag, double *dmax, double *dmin, double *dtrue, char *serr);

ext_def( double ) swe_deltat(double tjd);
ext_def(double) swe_deltat_ex(double tjd, int32 iflag, char *serr);

ext_def(int32) swe_time_equ(double tjd, double *te, char *serr);
ext_def(int32) swe_lmt_to_lat(double tjd_lmt, double geolon, double *tjd_lat, char *serr);
ext_def(int32) swe_lat_to_lmt(double tjd_lat, double geolon, double *tjd_lmt, char *serr);

ext_def( double ) swe_sidtime0(double tjd_ut, double eps, double nut);
ext_def( double ) swe_sidtime(double tjd_ut);
ext_def( void ) swe_set_interpolate_nut(AS_BOOL do_interpolate);

ext_def( void ) swe_cotrans(double *xpo, double *xpn, double eps);
ext_def( void ) swe_cotrans_sp(double *xpo, double *xpn, double eps);

ext_def( double ) swe_get_tid_acc(void);
ext_def( void ) swe_set_tid_acc(double t_acc);

ext_def (void) swe_set_delta_t_userdef(double dt);

ext_def( double ) swe_degnorm(double x);
ext_def( double ) swe_radnorm(double x);
ext_def( double ) swe_rad_midp(double x1, double x0);
ext_def( double ) swe_deg_midp(double x1, double x0);

ext_def( void ) swe_split_deg(double ddeg, int32 roundflag, int32 *ideg, int32 *imin, int32 *isec, double *dsecfr, int32 *isgn);

ext_def( centisec ) swe_csnorm(centisec p);

ext_def( centisec ) swe_difcsn (centisec p1, centisec p2);

ext_def( double ) swe_difdegn (double p1, double p2);

ext_def( centisec ) swe_difcs2n(centisec p1, centisec p2);

ext_def( double ) swe_difdeg2n(double p1, double p2);
ext_def( double ) swe_difrad2n(double p1, double p2);

ext_def( centisec ) swe_csroundsec(centisec x);

ext_def( int32 ) swe_d2l(double x);

ext_def( int ) swe_day_of_week(double jd);

ext_def( char *) swe_cs2timestr(CSEC t, int sep, AS_BOOL suppressZero, char *a);

ext_def( char *) swe_cs2lonlatstr(CSEC t, char pchar, char mchar, char *s);

ext_def( char *) swe_cs2degstr(CSEC t, char *a);

#endif

#endif

#ifdef __cplusplus
}
#endif
