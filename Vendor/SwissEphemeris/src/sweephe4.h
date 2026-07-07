# ifndef _EPHE_INCLUDED
# define _EPHE_INCLUDED
#ifdef __cplusplus
extern "C" {
#endif

#define PLACALC_SUN	0
#define PLACALC_EARTH	0
#define PLACALC_MOON	1
#define PLACALC_MERCURY 2
#define PLACALC_VENUS	3
#define PLACALC_MARS	4
#define PLACALC_JUPITER	5
#define PLACALC_SATURN	6
#define PLACALC_URANUS	7
#define PLACALC_NEPTUNE 8
#define PLACALC_PLUTO	9
#define PLACALC_LASTPLANET PLUTO
#define PLACALC_MEAN_NODE  10
#define PLACALC_TRUE_NODE  11
#define PLACALC_CHIRON	   12
#define PLACALC_LILITH	   13
#define PLACALC_CALC_N	  14
#define PLACALC_CERES     14
#define PLACALC_PALLAS    15
#define PLACALC_JUNO      16
#define PLACALC_VESTA     17
#define PLACALC_EARTHHEL  18
#define PLACALC_PFORTUNAE 19

# define PLACALC_AC	   20
# define PLACALC_ASC	   20
# define PLACALC_MC	   21
# define PLACALC_CALC_N_MC  22

# define EP4_BLOCKSIZE  sizeof(struct ep4)
# if HPUNIX
#   define EP4_PATH "/home/ephe/"
# else
#   define EP4_PATH  "ephe\\"
# endif
# define EP4_FILE	"sep4_"
# define EP4_NDAYS	10000L
# define NDB		10L

# define EP_NP (PLACALC_CHIRON + 3)

# define EP_ALL_PLANETS  ((1 << (PLACALC_CHIRON + 1)) - 1)
# define EP_CALC_N	(PLACALC_CHIRON+1)
# define EP_ECL_INDEX  (PLACALC_CHIRON + 1)
# define EP_NUT_INDEX  (PLACALC_CHIRON + 2)
# define EP_ECL_BIT   (1 << EP_ECL_INDEX)
# define EP_NUT_BIT   (1 << EP_NUT_INDEX)
# define EP_ALL_BITS   (EP_ALL_PLANETS|EP_ECL_BIT|EP_NUT_BIT)

# define EP_BIT_SPEED	16
# define EP_BIT_MUST_USE_EPHE 256

struct elon  {
	short	p0m;
	short	p0s;
	short	pd1m;
	short	pd1s;
	short	pd2[NDB-2];

	};

struct ep4  {
  short   j_10000;
  short   j_rest;

  short	ecl0m;
  short ecl0s;
  short ecld1[NDB-1];
  short	nuts[NDB];
  struct  elon elo[PLACALC_CHIRON +1];
};

extern FILE *ephfp;

extern centisec *ephread(double jd, int plalist, int flag, char *errtext);

extern double *dephread2(double jd, int plalist, int flag, char *errtext);

extern int eph4_posit (int jlong, AS_BOOL writeflag, char *errtext);

extern int ephe_plac2swe(int p);

extern void shortreorder (UCHAR *p, int n);

#ifdef __cplusplus
}
#endif

# endif
