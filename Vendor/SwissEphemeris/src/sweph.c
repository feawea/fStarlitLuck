#include <string.h>
#include <ctype.h>
#if MSDOS
#include <tchar.h>
#include <windows.h>
#endif
#include "swejpl.h"
#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"

#ifdef _MSC_VER
#define CMP_CALL_CONV __cdecl
#else
#define CMP_CALL_CONV
#endif

#define IS_PLANET 		0
#define IS_MOON			1
#define IS_ANY_BODY		2
#define IS_MAIN_ASTEROID	3

#define DO_SAVE			TRUE
#define NO_SAVE			FALSE

#define SEFLG_EPHMASK	(SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)
#define SEFLG_COORDSYS  (SEFLG_EQUATORIAL | SEFLG_XYZ | SEFLG_RADIANS)

struct meff_ele {double r,m;};

TLS struct swe_data swed = {FALSE,
                            FALSE,
                            NULL,
			    "",
			    "",
			    0,
			    0,
			    FALSE,
			    FALSE,
			    FALSE,
			    0.0, 0.0, 0.0, 0.0,
			    0,
			    0.0,
			    FALSE,
			    FALSE,
			    FALSE,
			    FALSE,
			    0.0,
			    0.0,
			    0.0,
			    0.0,
			    "",
			    0,
			    "",
			    NULL,
			    NULL,
			    0,
			    {0,0,0,0,0,0,0,0,},
			    };

static const char *ayanamsa_name[] = {
   "Fagan/Bradley",
   "Lahiri",
   "De Luce",
   "Raman",
   "Usha/Shashi",
   "Krishnamurti",
   "Djwhal Khul",
   "Yukteshwar",
   "J.N. Bhasin",
   "Babylonian/Kugler 1",
   "Babylonian/Kugler 2",
   "Babylonian/Kugler 3",
   "Babylonian/Huber",
   "Babylonian/Eta Piscium",
   "Babylonian/Aldebaran = 15 Tau",
   "Hipparchos",
   "Sassanian",
   "Galact. Center = 0 Sag",
   "J2000",
   "J1900",
   "B1950",
   "Suryasiddhanta",
   "Suryasiddhanta, mean Sun",
   "Aryabhata",
   "Aryabhata, mean Sun",
   "SS Revati",
   "SS Citra",
   "True Citra",
   "True Revati",
   "True Pushya (PVRN Rao)",
   "Galactic Center (Gil Brand)",
   "Galactic Equator (IAU1958)",
   "Galactic Equator",
   "Galactic Equator mid-Mula",
   "Skydram (Mardyks)",
   "True Mula (Chandra Hari)",
   "Dhruva/Gal.Center/Mula (Wilhelm)",
   "Aryabhata 522",
   "Babylonian/Britton",
   "\"Vedic\"/Sheoran",
   "Cochrane (Gal.Center = 0 Cap)",
   "Galactic Equator (Fiorenza)",
   "Vettius Valens",
   "Lahiri 1940",
   "Lahiri VP285",
   "Krishnamurti-Senthilathiban",
   "Lahiri ICRC",

};
static const int pnoint2jpl[]   = PNOINT2JPL;

static const int pnoext2int[] = {SEI_SUN, SEI_MOON, SEI_MERCURY, SEI_VENUS, SEI_MARS, SEI_JUPITER, SEI_SATURN, SEI_URANUS, SEI_NEPTUNE, SEI_PLUTO, 0, 0, 0, 0, SEI_EARTH, SEI_CHIRON, SEI_PHOLUS, SEI_CERES, SEI_PALLAS, SEI_JUNO, SEI_VESTA, };

static int32 swecalc(double tjd, int ipl, int iplmoon, int32 iflag, double *x, char *serr);
static int do_fread(void *targ, int size, int count, int corrsize,
		    FILE *fp, int32 fpos, int freord, int fendian, int ifno,
		    char *serr);
static int get_new_segment(double tjd, int ipli, int ifno, char *serr);
static int main_planet(double tjd, int ipli, int iplmoon, int32 epheflag, int32 iflag,
		       char *serr);
static int main_planet_bary(double tjd, int ipli, int32 epheflag, int32 iflag,
		AS_BOOL do_save,
		double *xp, double *xe, double *xs, double *xm,
		char *serr);
static int sweplan(double tjd, int ipli, int ifno, int32 iflag, AS_BOOL do_save,
		   double *xp, double *xpe, double *xps, double *xpm,
		   char *serr);
static int swemoon(double tjd, int32 iflag, AS_BOOL do_save, double *xp, char *serr);
static int sweph(double tjd, int ipli, int ifno, int32 iflag, double *xsunb, AS_BOOL do_save,
		double *xp, char *serr);
static int jplplan(double tjd, int ipli, int32 iflag, AS_BOOL do_save,
		   double *xp, double *xpe, double *xps, char *serr);
static void rot_back(int ipl);
static int read_const(int ifno, char *serr);
static void embofs(double *xemb, double *xmoon);
static int app_pos_etc_plan(int ipli, int iplmoon, int32 iflag, char *serr);
static int app_pos_etc_plan_osc(int ipl, int ipli, int32 iflag, char *serr);
static int app_pos_etc_sun(int32 iflag, char *serr);
static int app_pos_etc_moon(int32 iflag, char *serr);
static int app_pos_etc_sbar(int32 iflag, char *serr);
extern int swi_plan_for_osc_elem(int32 iflag, double tjd, double *xx);
static void swi_close_keep_topo_etc(void);
static int app_pos_etc_mean(int ipl, int32 iflag, char *serr);
static void nut_matrix(struct nut *nu, struct epsilon *oec);
static void calc_epsilon(double tjd, int32 iflag, struct epsilon *e);
static int lunar_osc_elem(double tjd, int ipl, int32 iflag, char *serr);
static int intp_apsides(double tjd, int ipl, int32 iflag, char *serr);
static double meff(double r);
static void denormalize_positions(double *x0, double *x1, double *x2);
static void calc_speed(double *x0, double *x1, double *x2, double dt);
static int32 plaus_iflag(int32 iflag, int32 ipl, double tjd, char *serr);
static int app_pos_rest(struct plan_data *pdp, int32 iflag,
    double *xx, double *x2000, struct epsilon *oe, char *serr);
static int open_jpl_file(double *ss, char *fname, char *fpath, char *serr);
static void free_planets(void);

#ifdef TRACE
static void trace_swe_calc(int param, double tjd, int ipl, int32 iflag, double *xx, char *serr);
static void trace_swe_fixstar(int swtch, char *star, double tjd, int32 iflag, double *xx, char *serr);
static void trace_swe_get_planet_name(int swtch, int ipl, char *s);
#endif

char *CALL_CONV swe_version(char *s)
{
  strcpy(s, SE_VERSION);
  return s;
}

#ifndef NO_SWE_GLP
#if MSDOS
HANDLE dllhandle = NULL;

#else
#ifdef __GNUC__

#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <dlfcn.h>
  static Dl_info dli;
#endif
#endif

char *CALL_CONV swe_get_library_path(char *s)
{
  size_t bytes;
  size_t len;
  *s = '\0';
#if !defined(__APPLE)
  len = AS_MAXCH;
  bytes = 0;
#if MSDOS
  bytes = GetModuleFileName((HMODULE) dllhandle, (TCHAR*) s, (DWORD) len);
#else
  #ifdef __GNUC__
    if (dladdr((void *)swe_version, &dli) != 0) {
      if (strlen(dli.dli_fname) >= len) {
	strncpy(s, dli.dli_fname, len);
	s[len] = '\0';
      } else{
	strcpy(s, dli.dli_fname);
      }
      bytes = strlen(s);
    } else {
      bytes = readlink("/proc/self/exe", s, len);
    }
  #else
    bytes = readlink("/proc/self/exe", s, len);
  #endif
#endif
  s[bytes] = '\0';
#endif
  return s;
}
#else

char *CALL_CONV swe_get_library_path(char *s)
{
  *s = '\0';
  return s;
}
#endif

int32 CALL_CONV swe_calc(double tjd, int ipl, int32 iflag,
	double *xx, char *serr)
{
  int i, j;
  int32 iplmoon = 0, iflgsave = iflag;
  int32 epheflag;
  AS_BOOL use_speed3 = FALSE;
  struct save_positions *sd;
  double x[6], *xs, x0[24], x2[24];
  double dt;
  if (serr != NULL)
    *serr = '\0';
#ifdef TRACE
#ifdef FORCE_IFLAG

  static TLS int force_flag = 0;
  static TLS int32 iflag_forced = 0;
  static TLS int force_flag_checked = 0;
  FILE *fp;
  char s[AS_MAXCH], *sp;
  memset(x, 0, sizeof(double) * 6);

  if (!force_flag_checked) {
    if ((fp = fopen(fname_force_flg, BFILE_R_ACCESS)) != NULL) {
      force_flag = 1;
      fgets(s, AS_MAXCH, fp);
      if ((sp = strchr(s, '\n')) != NULL)
	*sp = '\0';
      iflag_forced = atol(s);
      fclose(fp);
    }
    force_flag_checked = 1;
  }
  if (force_flag)
    iflag |= iflag_forced;
#endif
  swi_open_trace(serr);
  trace_swe_calc(1, tjd, ipl, iflag, xx, NULL);
#endif

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;

  epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag & SEFLG_MOSEPH) {
    epheflag = SEFLG_MOSEPH;
  } else if (epheflag & SEFLG_JPLEPH) {
    epheflag = SEFLG_JPLEPH;
  } else  {
    epheflag = SEFLG_SWIEPH;
  }
  if (swi_init_swed_if_start() == 1 && !(epheflag & SEFLG_MOSEPH) && serr != NULL) {
    strcpy(serr, "Please call swe_set_ephe_path() or swe_set_jplfile() before calling swe_calc() or swe_calc_ut()");
  }
  if (swed.last_epheflag != epheflag) {
    free_planets();

    if (ipl != SE_ECL_NUT) {
      if (swed.jpl_file_is_open) {
	swi_close_jpl_file();
	swed.jpl_file_is_open = FALSE;
      }
      for (i = 0; i < SEI_NEPHFILES; i ++) {
	if (swed.fidat[i].fptr != NULL)
	  fclose(swed.fidat[i].fptr);
	memset((void *) &swed.fidat[i], 0, sizeof(struct file_data));
      }
      swed.last_epheflag = epheflag;
    }
  }

  if ((iflag & SEFLG_SPEED3) && (iflag & SEFLG_SPEED))
    iflag = iflag & ~SEFLG_SPEED3;
  if (iflag & SEFLG_SPEED3)
    use_speed3 = TRUE;

  if ((iflag & SEFLG_SPEED) && (iflag & SEFLG_TOPOCTR) && !(iflag & SEFLG_NOABERR))
    use_speed3 = TRUE;

  if ((iflag & SEFLG_XYZ) && (iflag & SEFLG_RADIANS))
    iflag = iflag & ~SEFLG_RADIANS;

  if ((iflag & SEFLG_CENTER_BODY) && ipl <= SE_PLUTO && (iflag & SEFLG_TEST_PLMOON) != SEFLG_TEST_PLMOON) {
    iplmoon = ipl * 100 + 9099;
  }

  if (ipl >= SE_PLMOON_OFFSET && ipl < SE_AST_OFFSET && (iflag & SEFLG_TEST_PLMOON) != SEFLG_TEST_PLMOON) {
    iplmoon = ipl;
    ipl = (int) ((ipl - 9000) / 100);
    iflag |= SEFLG_CENTER_BODY;
  }

  if ((iflag & SEFLG_CENTER_BODY) && ipl <= SE_MARS && (iplmoon % 100) == 99) {
    iplmoon = 0;
    iflag &= ~SEFLG_CENTER_BODY;
  }
  if ((iflag & SEFLG_CENTER_BODY) || iplmoon > 0)
    swi_force_app_pos_etc();

  if (ipl < SE_NPLANETS && ipl >= SE_SUN) {
    sd = &swed.savedat[ipl];

  } else {

    sd = &swed.savedat[SE_NPLANETS];
  }

  if (sd->tsave == tjd && tjd != 0 && ipl == sd->ipl && iplmoon == 0) {
    if ((sd->iflgsave & ~SEFLG_COORDSYS) == (iflag & ~SEFLG_COORDSYS))
      goto end_swe_calc;
  }

  if (!use_speed3) {

    sd->tsave = tjd;
    sd->ipl = ipl;
    if ((sd->iflgsave = swecalc(tjd, ipl, iplmoon, iflag, sd->xsaves, serr)) == ERR)
      goto return_error;
  } else {

    sd->tsave = tjd;
    sd->ipl = ipl;
    switch(ipl) {
      case SE_MOON:
	dt = MOON_SPEED_INTV;
	break;
      case SE_OSCU_APOG:
      case SE_TRUE_NODE:

	dt = NODE_CALC_INTV_MOSH;
	break;
      default:
	dt = PLAN_SPEED_INTV;
	break;
    }
    if ((sd->iflgsave = swecalc(tjd-dt, ipl, iplmoon, iflag, x0, serr)) == ERR)
      goto return_error;
    if ((sd->iflgsave = swecalc(tjd+dt, ipl, iplmoon, iflag, x2, serr)) == ERR)
      goto return_error;
    if ((sd->iflgsave = swecalc(tjd, ipl, iplmoon, iflag, sd->xsaves, serr)) == ERR)
      goto return_error;
    denormalize_positions(x0, sd->xsaves, x2);
    calc_speed(x0, sd->xsaves, x2, dt);
  }
  end_swe_calc:
  if (iflag & SEFLG_EQUATORIAL) {
    xs = sd->xsaves+12;
  } else {
    xs = sd->xsaves;
  }
  if (iflag & SEFLG_XYZ)
    xs = xs+6;
  if (ipl == SE_ECL_NUT) {
    i = 4;
  } else {
    i = 3;
  }
  for (j = 0; j < i; j++)
    x[j] = *(xs + j);
  for (j = i; j < 6; j++)
    x[j] = 0;
  if (iflag & (SEFLG_SPEED3 | SEFLG_SPEED)) {
    for (j = 3; j < 6; j++)
      x[j] = *(xs + j);
  }
#if 1
  if (iflag & SEFLG_RADIANS) {
    if (ipl == SE_ECL_NUT) {
      for (j = 0; j < 4; j++)
        x[j] *= DEGTORAD;
    } else {
      for (j = 0; j < 2; j++)
        x[j] *= DEGTORAD;
      if (iflag & (SEFLG_SPEED3 | SEFLG_SPEED)) {
        for (j = 3; j < 5; j++)
	  x[j] *= DEGTORAD;
      }
    }
  }
#endif
  for (i = 0; i <= 5; i++)
    xx[i] = x[i];

  iflag = sd->iflgsave & ~SEFLG_COORDSYS;

  iflag |= (iflgsave & SEFLG_COORDSYS);

  if ((iflgsave & SEFLG_EPHMASK) == 0)
    iflag = iflag & ~SEFLG_DEFAULTEPH;
#ifdef TRACE
  trace_swe_calc(2, tjd, ipl, iflag, xx, serr);
#endif
  return iflag;
return_error:
  for (i = 0; i <= 5; i++)
    xx[i] = 0;
#ifdef TRACE
  trace_swe_calc(2, tjd, ipl, iflag, xx, serr);
#endif
  return ERR;
}

int32 CALL_CONV swe_calc_ut(double tjd_ut, int32 ipl, int32 iflag,
	double *xx, char *serr)
{
  double deltat;
  int32 retval = OK;
  int32 epheflag = 0;
  iflag = plaus_iflag(iflag, ipl, tjd_ut, serr);
  epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag == 0) {
    epheflag = SEFLG_SWIEPH;
    iflag |= SEFLG_SWIEPH;
  }
  deltat = swe_deltat_ex(tjd_ut, iflag, serr);
  retval = swe_calc(tjd_ut + deltat, ipl, iflag, xx, serr);

  if ((retval & SEFLG_EPHMASK) != epheflag) {
    deltat = swe_deltat_ex(tjd_ut, retval, NULL);
    retval = swe_calc(tjd_ut + deltat, ipl, iflag, xx, NULL);
  }
  return retval;
}

static int32 swecalc(double tjd, int ipl, int32 iplmoon, int32 iflag, double *x, char *serr)
{
  int i;
  int ipli, ipli_ast, ifno;
  int retc;
  int32 epheflag = SEFLG_DEFAULTEPH;
  struct plan_data *pdp;
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  struct plan_data *ndp;
  double *xp, *xp2;
  double ss[3];
  char serr2[AS_MAXCH];

  serr2[0] = '\0';

  iflag = plaus_iflag(iflag, ipl, tjd, serr);

  if (iflag & SEFLG_MOSEPH)
    epheflag = SEFLG_MOSEPH;
  if (iflag & SEFLG_SWIEPH)
    epheflag = SEFLG_SWIEPH;
  if (iflag & SEFLG_JPLEPH)
    epheflag = SEFLG_JPLEPH;

  if ((iflag & SEFLG_BARYCTR) && (iflag & SEFLG_MOSEPH)) {
    if (serr != NULL)
      strcpy(serr, "barycentric Moshier positions are not supported.");
    return ERR;
  }
  if (epheflag != SEFLG_MOSEPH && !swed.ephe_path_is_set && !swed.jpl_file_is_open)
    swe_set_ephe_path(NULL);
  if ((iflag & SEFLG_SIDEREAL) && !swed.ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);

  swi_check_ecliptic(tjd, iflag);

  swi_check_nutation(tjd, iflag);

  if (ipl == SE_ECL_NUT) {
    x[0] = swed.oec.eps + swed.nut.nutlo[1];
    x[1] = swed.oec.eps;
    x[2] = swed.nut.nutlo[0];
    x[3] = swed.nut.nutlo[1];

      for (i = 0; i <= 3; i++)
	x[i] *= RADTODEG;
    return(iflag);

  } else if (ipl == SE_MOON) {

    ipli = SEI_MOON;
    pdp = &swed.pldat[ipli];
    xp = pdp->xreturn;
    switch(epheflag) {
      case SEFLG_JPLEPH:
	retc = jplplan(tjd, ipli, iflag, DO_SAVE, NULL, NULL, NULL, serr);

	if (retc == ERR)
	  goto return_error;

        if (retc == NOT_AVAILABLE) {
	  iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_SWIEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \ntrying Swiss Eph; ");
	  goto sweph_moon;
	} else if (retc == BEYOND_EPH_LIMITS) {
	  if (tjd > MOSHLUEPH_START && tjd < MOSHLUEPH_END) {
	    iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_MOSEPH;
	    if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	      strcat(serr, " \nusing Moshier Eph; ");
	    goto moshier_moon;
	  } else
	    goto return_error;
	}
	break;
      case SEFLG_SWIEPH:
	sweph_moon:
	retc = sweplan(tjd, ipli, SEI_FILE_MOON, iflag, DO_SAVE,
			NULL, NULL, NULL, NULL, serr);
	if (retc == ERR)
	  goto return_error;

        if (retc == NOT_AVAILABLE) {
	  if (tjd > MOSHLUEPH_START && tjd < MOSHLUEPH_END) {
	    iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
	    if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	      strcat(serr, " \nusing Moshier eph.; ");
	    goto moshier_moon;
	  } else
	    goto return_error;
	}
	break;
      case SEFLG_MOSEPH:
	moshier_moon:
        retc = swi_moshmoon(tjd, DO_SAVE, NULL, serr);
	if (retc == ERR)
	  goto return_error;

	retc = swi_moshplan(tjd, SEI_EARTH, DO_SAVE, NULL, NULL, serr);
	if (retc == ERR)
	  goto return_error;
	break;
      default:
	break;
    }

    if ((retc = app_pos_etc_moon(iflag, serr)) != OK)
      goto return_error;

  } else if (ipl == SE_SUN && (iflag & SEFLG_BARYCTR)) {

    ipli = SEI_SUN;
    xp = pedp->xreturn;
    switch(epheflag) {
      case SEFLG_JPLEPH:

	if (!swed.jpl_file_is_open) {
	  retc = open_jpl_file(ss, swed.jplfnam, swed.ephepath, serr);
	  if (retc != OK)
	    goto sweph_sbar;
	}
	retc = swi_pleph(tjd, J_SUN, J_SBARY, psdp->x, serr);
	if (retc == ERR || retc == BEYOND_EPH_LIMITS) {
	  swi_close_jpl_file();
	  swed.jpl_file_is_open = FALSE;
	  goto return_error;
	}

        if (retc == NOT_AVAILABLE) {
	  iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_SWIEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \ntrying Swiss Eph; ");
	  goto sweph_sbar;
	}
	psdp->teval = tjd;
	break;
      case SEFLG_SWIEPH:
	sweph_sbar:

	retc = sweplan(tjd, SEI_EARTH, SEI_FILE_PLANET, iflag, DO_SAVE, NULL, NULL, NULL, NULL, serr);
#if 1
	if (retc == ERR || retc == NOT_AVAILABLE)
	  goto return_error;
#else

	if (retc == ERR)
	  goto return_error;

        if (retc == NOT_AVAILABLE) {
	  if (tjd > MOSHLUEPH_START && tjd < MOSHLUEPH_END) {
	    iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
	    if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	      strcat(serr, " \nusing Moshier; ");
	    goto moshier_sbar;
	  } else
	    goto return_error;
	}
#endif
	psdp->teval = tjd;

	break;
      default:
	return ERR;
	break;
    }

    if ((retc = app_pos_etc_sbar(iflag, serr)) != OK)
      goto return_error;

    iflag = pedp->xflgs;

    pedp->xflgs = -1;

  } else if (ipl == SE_SUN
	  || ipl == SE_MERCURY
	  || ipl == SE_VENUS
	  || ipl == SE_MARS
	  || ipl == SE_JUPITER
	  || ipl == SE_SATURN
	  || ipl == SE_URANUS
	  || ipl == SE_NEPTUNE
	  || ipl == SE_PLUTO
	  || ipl == SE_EARTH) {
    if (iflag & SEFLG_HELCTR) {
      if (ipl == SE_SUN) {

	for (i = 0; i < 24; i++)
	  x[i] = 0;
	return iflag;
      }
    } else if (iflag & SEFLG_BARYCTR) {
      ;
    } else {
      if (ipl == SE_EARTH) {

	for (i = 0; i < 24; i++)
	  x[i] = 0;
	return iflag;
      }
    }

    ipli = pnoext2int[ipl];
    pdp = &swed.pldat[ipli];
    xp = pdp->xreturn;
    retc = main_planet(tjd, ipli, iplmoon, epheflag, iflag, serr);
    if (retc == ERR)
      goto return_error;

    iflag = pdp->xflgs;

  } else if (ipl == SE_MEAN_NODE) {
    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {

      for (i = 0; i < 24; i++)
	x[i] = 0;
      return iflag;
    }
    ndp = &swed.nddat[SEI_MEAN_NODE];
    xp = ndp->xreturn;
    xp2 = ndp->x;
    retc = swi_mean_node(tjd, xp2, serr);
    if (retc == ERR)
      goto return_error;

    retc = swi_mean_node(tjd - MEAN_NODE_SPEED_INTV, xp2+3, serr);
    if (retc == ERR)
      goto return_error;
    xp2[3] = swe_difrad2n(xp2[0], xp2[3]) / MEAN_NODE_SPEED_INTV;
    xp2[4] = xp2[5] = 0;
    ndp->teval = tjd;
    ndp->xflgs = -1;

    if ((retc = app_pos_etc_mean(SEI_MEAN_NODE, iflag, serr)) != OK)
      goto return_error;

    if (!(iflag & SEFLG_SIDEREAL) && !(iflag & SEFLG_J2000)) {
      ndp->xreturn[1] = 0.0;
      ndp->xreturn[4] = 0.0;
      ndp->xreturn[5] = 0.0;
      ndp->xreturn[8] = 0.0;
      ndp->xreturn[11] = 0.0;
    }
    if (retc == ERR)
      goto return_error;

  } else if (ipl == SE_MEAN_APOG) {
    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {

      for (i = 0; i < 24; i++)
	x[i] = 0;
      return iflag;
    }
    ndp = &swed.nddat[SEI_MEAN_APOG];
    xp = ndp->xreturn;
    xp2 = ndp->x;
    retc = swi_mean_apog(tjd, xp2, serr);
    if (retc == ERR)
      goto return_error;

    retc = swi_mean_apog(tjd - MEAN_NODE_SPEED_INTV, xp2+3, serr);
    if (retc == ERR)
      goto return_error;
    for(i = 0; i <= 1; i++)
      xp2[3+i] = swe_difrad2n(xp2[i], xp2[3+i]) / MEAN_NODE_SPEED_INTV;
    xp2[5] = 0;
    ndp->teval = tjd;
    ndp->xflgs = -1;

    if ((retc = app_pos_etc_mean(SEI_MEAN_APOG, iflag, serr)) != OK)
      goto return_error;

    ndp->xreturn[5] = 0.0;
    if (retc == ERR)
      goto return_error;

  } else if (ipl == SE_TRUE_NODE) {
    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {

      for (i = 0; i < 24; i++)
	x[i] = 0;
      return iflag;
    }
    ndp = &swed.nddat[SEI_TRUE_NODE];
    xp = ndp->xreturn;
    retc = lunar_osc_elem(tjd, SEI_TRUE_NODE, iflag, serr);
    iflag = ndp->xflgs;

    if (!(iflag & SEFLG_SIDEREAL) && !(iflag & SEFLG_J2000)) {
      ndp->xreturn[1] = 0.0;
      ndp->xreturn[4] = 0.0;
      ndp->xreturn[8] = 0.0;
      ndp->xreturn[11] = 0.0;
    }
    if (retc == ERR)
      goto return_error;

  } else if (ipl == SE_OSCU_APOG) {
    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {

      for (i = 0; i < 24; i++)
	x[i] = 0;
      return iflag;
    }
    ndp = &swed.nddat[SEI_OSCU_APOG];
    xp = ndp->xreturn;
    retc = lunar_osc_elem(tjd, SEI_OSCU_APOG, iflag, serr);
    iflag = ndp->xflgs;
    if (retc == ERR)
      goto return_error;

  } else if (ipl == SE_INTP_APOG) {
    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {

      for (i = 0; i < 24; i++)
	x[i] = 0;
      return iflag;
    }
    if (tjd < MOSHLUEPH_START || tjd > MOSHLUEPH_END) {
      for (i = 0; i < 24; i++)
	x[i] = 0;
      if (serr != NULL)
	sprintf(serr, "Interpolated apsides are restricted to JD %8.1f - JD %8.1f",
		MOSHLUEPH_START, MOSHLUEPH_END);
      return ERR;
    }
    ndp = &swed.nddat[SEI_INTP_APOG];
    xp = ndp->xreturn;
    retc = intp_apsides(tjd, SEI_INTP_APOG, iflag, serr);
    iflag = ndp->xflgs;
    if (retc == ERR)
      goto return_error;

  } else if (ipl == SE_INTP_PERG) {
    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {

      for (i = 0; i < 24; i++)
	x[i] = 0;
      return iflag;
    }
    if (tjd < MOSHLUEPH_START || tjd > MOSHLUEPH_END) {
      for (i = 0; i < 24; i++)
	x[i] = 0;
      if (serr != NULL)
	sprintf(serr, "Interpolated apsides are restricted to JD %8.1f - JD %8.1f",
		MOSHLUEPH_START, MOSHLUEPH_END);
      return ERR;
    }
    ndp = &swed.nddat[SEI_INTP_PERG];
    xp = ndp->xreturn;
    retc = intp_apsides(tjd, SEI_INTP_PERG, iflag, serr);
    iflag = ndp->xflgs;
    if (retc == ERR)
      goto return_error;

  } else if (ipl == SE_CHIRON
    || ipl == SE_PHOLUS
    || ipl == SE_CERES
    || ipl == SE_PALLAS
    || ipl == SE_JUNO
    || ipl == SE_VESTA
    || ipl > SE_PLMOON_OFFSET
    || ipl > SE_AST_OFFSET
    ) {

    if (ipl < SE_NPLANETS) {
      ipli = pnoext2int[ipl];
    } else if (ipl <= SE_AST_OFFSET + MPC_VESTA && ipl > SE_AST_OFFSET) {
      ipli = SEI_CERES + ipl - SE_AST_OFFSET - 1;
      ipl = SE_CERES + ipl - SE_AST_OFFSET - 1;
    } else {
      ipli = SEI_ANYBODY;
    }
    if (ipli == SEI_ANYBODY) {
      ipli_ast = ipl;
    } else {
      ipli_ast = ipli;
    }
    pdp = &swed.pldat[ipli];
    xp = pdp->xreturn;
    if (ipli_ast > SE_AST_OFFSET) {
      ifno = SEI_FILE_ANY_AST;
    } else if (ipli_ast > SE_PLMOON_OFFSET) {
      ifno = SEI_FILE_ANY_AST;
    } else {
      ifno = SEI_FILE_MAIN_AST;
    }
    if (ipli == SEI_CHIRON && (tjd < CHIRON_START || tjd > CHIRON_END)) {
      if (serr != NULL)
	sprintf(serr, "Chiron's ephemeris is restricted to JD %8.1f - JD %8.1f",
		CHIRON_START, CHIRON_END);
      return ERR;
    }
    if (ipli == SEI_PHOLUS && (tjd < PHOLUS_START || tjd > PHOLUS_END)) {
      if (serr != NULL)
	sprintf(serr,
		"Pholus's ephemeris is restricted to JD %8.1f - JD %8.1f",
		PHOLUS_START, PHOLUS_END);
      return ERR;
    }
  do_asteroid:

    retc = main_planet(tjd, SEI_EARTH, 0, epheflag, iflag, serr);
    if (retc == ERR)
      goto return_error;

    iflag = swed.pldat[SEI_EARTH].xflgs;

    if (serr != NULL) {
      strcpy(serr2, serr);
      *serr = '\0';
    }

    retc = sweph(tjd, ipli_ast, ifno, iflag, psdp->x, DO_SAVE, NULL, serr);
    if (retc == ERR || retc == NOT_AVAILABLE)
      goto return_error;
    retc = app_pos_etc_plan(ipli_ast, 0, iflag, serr);
    if (retc == ERR)
      goto return_error;

    if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS) {
      if (epheflag != SEFLG_MOSEPH) {
	iflag = (iflag & ~SEFLG_EPHMASK) | SEFLG_MOSEPH;
	epheflag = SEFLG_MOSEPH;
	if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	  strcat(serr, "\nusing Moshier eph.; ");
	goto do_asteroid;
      } else
	goto return_error;
    }

    if (serr != NULL && *serr == '\0' && *serr2 != '\0') {
      strcpy(serr, "sun: ");
      serr2[AS_MAXCH-5] = '\0';
      strcat(serr, serr2);
    }

  } else if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX) {

    ipli = SEI_ANYBODY;
    pdp = &swed.pldat[ipli];
    xp = pdp->xreturn;
  do_fict_plan:

    retc = main_planet(tjd, SEI_EARTH, 0, epheflag, iflag, serr);

    iflag = swed.pldat[SEI_EARTH].xflgs;

    if (swi_osc_el_plan(tjd, pdp->x, ipl-SE_FICT_OFFSET, ipli, pedp->x, psdp->x, serr) != OK)
      goto return_error;
    if (retc == ERR)
      goto return_error;
    retc = app_pos_etc_plan_osc(ipl, ipli, iflag, serr);
    if (retc == ERR)
      goto return_error;

    if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS) {
      if (epheflag != SEFLG_MOSEPH) {
	iflag = (iflag & ~SEFLG_EPHMASK) | SEFLG_MOSEPH;
	epheflag = SEFLG_MOSEPH;
	if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	  strcat(serr, "\nusing Moshier eph.; ");
	goto do_fict_plan;
      } else
	goto return_error;
    }

  } else {
    if (serr != NULL) {
      sprintf(serr, "illegal planet number %d.", ipl);
    }
    goto return_error;
  }
  for (i = 0; i < 24; i++)
    x[i] = xp[i];
  return(iflag);

  return_error:;
  for (i = 0; i < 24; i++)
    x[i] = 0;
  return ERR;
}

static void free_planets(void)
{
  int i;

  for (i = 0; i < SEI_NPLANETS; i++) {
    if (swed.pldat[i].segp != NULL) {
      free((void *) swed.pldat[i].segp);
    }
    if (swed.pldat[i].refep != NULL) {
      free((void *) swed.pldat[i].refep);
    }
    memset((void *) &swed.pldat[i], 0, sizeof(struct plan_data));
  }
  for (i = 0; i <= SE_NPLANETS; i++)
    memset((void *) &swed.savedat[i], 0, sizeof(struct save_positions));

  for (i = 0; i < SEI_NNODE_ETC; i++) {
    memset((void *) &swed.nddat[i], 0, sizeof(struct plan_data));
  }
}

int32 swi_init_swed_if_start(void)
{

  if (!swed.swed_is_initialised) {
    memset((void *) &swed, 0, sizeof(struct swe_data));
    strcpy(swed.ephepath, SE_EPHE_PATH);
    strcpy(swed.jplfnam, SE_FNAME_DFT);
    swe_set_tid_acc(SE_TIDAL_AUTOMATIC);
    swed.swed_is_initialised = TRUE;
    return 1;
  }
  return 0;
}

static void swi_close_keep_topo_etc(void)
{
  int i;

  for (i = 0; i < SEI_NEPHFILES; i ++) {
    if (swed.fidat[i].fptr != NULL)
      fclose(swed.fidat[i].fptr);
    memset((void *) &swed.fidat[i], 0, sizeof(struct file_data));
  }
  free_planets();
  memset((void *) &swed.oec, 0, sizeof(struct epsilon));
  memset((void *) &swed.oec2000, 0, sizeof(struct epsilon));
  memset((void *) &swed.nut, 0, sizeof(struct nut));
  memset((void *) &swed.nut2000, 0, sizeof(struct nut));
  memset((void *) &swed.nutv, 0, sizeof(struct nut));
  memset((void *) &swed.astro_models, 0, SEI_NMODELS * sizeof(int32));

  swi_close_jpl_file();
  swed.jpl_file_is_open = FALSE;
  swed.jpldenum = 0;

  if (swed.fixfp != NULL) {
    fclose(swed.fixfp);
    swed.fixfp = NULL;
  }
  swe_set_tid_acc(SE_TIDAL_AUTOMATIC);
  swed.is_old_starfile = FALSE;
  swed.i_saved_planet_name = 0;
  *(swed.saved_planet_name) = '\0';
  swed.timeout = 0;
}

void CALL_CONV swe_close(void)
{
  int i;

  for (i = 0; i < SEI_NEPHFILES; i ++) {
    if (swed.fidat[i].fptr != NULL)
      fclose(swed.fidat[i].fptr);
    memset((void *) &swed.fidat[i], 0, sizeof(struct file_data));
  }
  free_planets();
  memset((void *) &swed.oec, 0, sizeof(struct epsilon));
  memset((void *) &swed.oec2000, 0, sizeof(struct epsilon));
  memset((void *) &swed.nut, 0, sizeof(struct nut));
  memset((void *) &swed.nut2000, 0, sizeof(struct nut));
  memset((void *) &swed.nutv, 0, sizeof(struct nut));
  memset((void *) &swed.astro_models, 0, SEI_NMODELS * sizeof(int32));

  swi_close_jpl_file();
  swed.jpl_file_is_open = FALSE;
  swed.jpldenum = 0;

  if (swed.fixfp != NULL) {
    fclose(swed.fixfp);
    swed.fixfp = NULL;
  }
  swe_set_tid_acc(SE_TIDAL_AUTOMATIC);
  swed.geopos_is_set = FALSE;
  swed.ayana_is_set = FALSE;
  swed.is_old_starfile = FALSE;
  swed.i_saved_planet_name = 0;
  *(swed.saved_planet_name) = '\0';
  memset((void *) &swed.topd, 0, sizeof(struct topo_data));
  memset((void *) &swed.sidd, 0, sizeof(struct sid_data));
  swed.timeout = 0;
  swed.last_epheflag = 0;
  if (swed.dpsi != NULL) {
    free(swed.dpsi);
    swed.dpsi = NULL;
  }
  if (swed.deps != NULL) {
    free(swed.deps);
    swed.deps = NULL;
  }
  if (swed.n_fixstars_records > 0) {
    free(swed.fixed_stars);
    swed.fixed_stars = NULL;
    swed.n_fixstars_real = 0;
    swed.n_fixstars_named = 0;
    swed.n_fixstars_records = 0;
  }

#ifdef TRACE
#define TRACE_CLOSE FALSE
  swi_open_trace(NULL);
  if (swi_fp_trace_c != NULL) {
    if (swi_trace_count < TRACE_COUNT_MAX) {
      fputs("\n[SWE_CLOSE]\n", swi_fp_trace_c);
      fputs("  swe_close();\n", swi_fp_trace_c);
#if TRACE_CLOSE
      fputs("}\n", swi_fp_trace_c);
#endif
      fflush(swi_fp_trace_c);
    }
#if TRACE_CLOSE
    fclose(swi_fp_trace_c);
#endif
  }
#if TRACE_CLOSE
  if (swi_fp_trace_out != NULL)
    fclose(swi_fp_trace_out);
  swi_fp_trace_c = NULL;
  swi_fp_trace_out = NULL;
#endif
#endif
}

void CALL_CONV swe_set_ephe_path(const char *path)
{
  int i, iflag;
  char s[AS_MAXCH];
  char *sp;
  double xx[6];

  swi_close_keep_topo_etc();
  swi_init_swed_if_start();
  swed.ephe_path_is_set = TRUE;

  if ((sp = getenv("SE_EPHE_PATH")) != NULL
    && strlen(sp) != 0
    && strlen(sp) <= AS_MAXCH-1-13) {
    strcpy(s, sp);
  } else if (path == NULL || *path == '\0') {
    strcpy(s, SE_EPHE_PATH);
  } else if (strlen(path) <= AS_MAXCH-1-13) {
    strcpy(s, path);
  } else {
    strcpy(s, SE_EPHE_PATH);
  }
  i = (int) strlen(s);
  if (*(s + i - 1) != *DIR_GLUE && *s != '\0')
    strcat(s, DIR_GLUE);
  strcpy(swed.ephepath, s);

  iflag = SEFLG_SWIEPH|SEFLG_J2000|SEFLG_TRUEPOS|SEFLG_ICRS;
  swed.last_epheflag = 2;
  swe_calc(J2000, SE_MOON, iflag, xx, NULL);
  if (swed.fidat[SEI_FILE_MOON].fptr != NULL) {
    swi_set_tid_acc(0, 0, swed.fidat[SEI_FILE_MOON].sweph_denum, NULL);
  }
#ifdef TRACE
  swi_open_trace(NULL);
  if (swi_trace_count < TRACE_COUNT_MAX) {
    if (swi_fp_trace_c != NULL) {
      fputs("\n[SWE_SET_EPHE_PATH]\n", swi_fp_trace_c);
      if (path == NULL) {
        fputs("  *s = '\\0';\n", swi_fp_trace_c);
      } else {
	fprintf(swi_fp_trace_c, "  strcpy(s, \"%s\");\n", path);
      }
      fputs("  swe_set_ephe_path(s);\n", swi_fp_trace_c);
      fputs("  printf(\"swe_set_ephe_path: path_in = \");", swi_fp_trace_c);
      fputs("  printf(s);\n", swi_fp_trace_c);
      fputs("  \tprintf(\"\\tpath_set = unknown to swetrace\\n\"); [ unknown to swetrace ]\n", swi_fp_trace_c);
      fflush(swi_fp_trace_c);
    }
    if (swi_fp_trace_out != NULL) {
      fputs("swe_set_ephe_path: path_in = ", swi_fp_trace_out);
      if (path != NULL)
	fputs(path, swi_fp_trace_out);
      fputs("\tpath_set = ", swi_fp_trace_out);
      fputs(s, swi_fp_trace_out);
      fputs("\n", swi_fp_trace_out);
      fflush(swi_fp_trace_out);
    }
  }
#endif
}

void load_dpsi_deps(void)
{
  FILE *fp;
  char s[AS_MAXCH];
  char *cpos[20];
  int n = 0, iyear, mjd = 0, mjdsv = 0;
  double dpsi, deps, TJDOFS = 2400000.5;
  if (swed.eop_dpsi_loaded > 0)
    return;
  fp = swi_fopen(-1, DPSI_DEPS_IAU1980_FILE_EOPC04, swed.ephepath, NULL);
  if (fp == NULL) {
    swed.eop_dpsi_loaded = ERR;
    return;
  }
  if ((swed.dpsi = (double *) calloc((size_t) SWE_DATA_DPSI_DEPS, sizeof(double))) == NULL) {
    swed.eop_dpsi_loaded = ERR;
    return;
  }
  if ((swed.deps = (double *) calloc((size_t) SWE_DATA_DPSI_DEPS, sizeof(double))) == NULL) {
    swed.eop_dpsi_loaded = ERR;
    return;
  }
  swed.eop_tjd_beg_horizons = DPSI_DEPS_IAU1980_TJD0_HORIZONS;
  while (fgets(s, AS_MAXCH, fp) != NULL) {
    swi_cutstr(s, " ", cpos, 16);
    if ((iyear = atoi(cpos[0])) == 0)
      continue;
    mjd = atoi(cpos[3]);

    if (mjdsv > 0 && mjd - mjdsv != 1) {

      swed.eop_dpsi_loaded = -2;
      fclose(fp);
      return;
    }
    if (n == 0)
      swed.eop_tjd_beg = mjd + TJDOFS;
    swed.dpsi[n] = atof(cpos[8]);
    swed.deps[n] = atof(cpos[9]);

    n++;
    mjdsv = mjd;
  }
  swed.eop_tjd_end = mjd + TJDOFS;
  swed.eop_dpsi_loaded = 1;
  fclose(fp);

  fp = swi_fopen(-1, DPSI_DEPS_IAU1980_FILE_FINALS, swed.ephepath, NULL);
  if (fp == NULL)
    return;
  while (fgets(s, AS_MAXCH, fp) != NULL) {
    mjd = atoi(s + 7);
    if (mjd + TJDOFS <= swed.eop_tjd_end)
      continue;
    if (n >= SWE_DATA_DPSI_DEPS)
      return;

    if (mjdsv > 0 && mjd - mjdsv != 1) {

      swed.eop_dpsi_loaded = -3;
      fclose(fp);
      return;
    }

    dpsi = atof(s + 168);
    deps = atof(s + 178);
    if (dpsi == 0) {

      dpsi = atof(s + 99);
      deps = atof(s + 118);
    }
    if (dpsi == 0) {
      swed.eop_dpsi_loaded = 2;

      fclose(fp);
      return;
    }
    swed.eop_tjd_end = mjd + TJDOFS;
    swed.dpsi[n] = dpsi / 1000.0;
    swed.deps[n] = deps / 1000.0;

    n++;
    mjdsv = mjd;
  }
  swed.eop_dpsi_loaded = 2;
  fclose(fp);
}

void CALL_CONV swe_set_jpl_file(const char *fname)
{
  char *sp, s[AS_MAXCH];
  int retc;
  double ss[3];

  swi_close_keep_topo_etc();
  swi_init_swed_if_start();

  if (strlen(fname) >= AS_MAXCH) {
     strncpy(s, fname, AS_MAXCH - 1);
     s[AS_MAXCH - 1] = '\0';
  } else {
    strcpy(s, fname);
  }
  sp = strrchr(s, (int) *DIR_GLUE);
  if (sp == NULL) {
    sp = s;
  } else {
    sp = sp + 1;
  }
  if (strlen(sp) >= AS_MAXCH)
    sp[AS_MAXCH - 1] = '\0';
  strcpy(swed.jplfnam, sp);

  retc = open_jpl_file(ss, swed.jplfnam, swed.ephepath, NULL);
  if (retc == OK) {
    if (swed.jpldenum >= 403) {

	load_dpsi_deps();
    }
  }
#ifdef TRACE
  swi_open_trace(NULL);
  if (swi_trace_count < TRACE_COUNT_MAX) {
    if (swi_fp_trace_c != NULL) {
      fputs("\n[SWE_SET_JPL_FILE]\n", swi_fp_trace_c);
      fprintf(swi_fp_trace_c, "  strcpy(s, \"%s\");\n", fname);
      fputs("  swe_set_jpl_file(s);\n", swi_fp_trace_c);
      fputs("  printf(\"swe_set_jpl_file: fname_in = \");", swi_fp_trace_c);
      fputs("  printf(s);\n", swi_fp_trace_c);
      fputs("  printf(\"\\tfname_set = unknown to swetrace\\n\");  [ unknown to swetrace ]\n", swi_fp_trace_c);
      fflush(swi_fp_trace_c);
    }
    if (swi_fp_trace_out != NULL) {
      fputs("swe_set_jpl_file: fname_in = ", swi_fp_trace_out);
      fputs(fname, swi_fp_trace_out);
      fputs("\tfname_set = ", swi_fp_trace_out);
      fputs(sp, swi_fp_trace_out);
      fputs("\n", swi_fp_trace_out);
      fflush(swi_fp_trace_out);
    }
  }
#endif
}

static void calc_epsilon(double tjd, int32 iflag, struct epsilon *e)
{
    e->teps = tjd;
    e->eps = swi_epsiln(tjd, iflag);
    e->seps = sin(e->eps);
    e->ceps = cos(e->eps);
}

static int main_planet(double tjd, int ipli, int iplmoon, int32 epheflag, int32 iflag,
		       char *serr)
{
  int retc;
  if ((iflag & SEFLG_CENTER_BODY)
    && ipli >= SE_MARS && ipli <= SE_PLUTO) {

    retc = sweph(tjd, iplmoon, SEI_FILE_ANY_AST, iflag, NULL, DO_SAVE, NULL, serr);
    if (retc == ERR || retc == NOT_AVAILABLE)
      return ERR;
  }
  switch(epheflag) {
    case SEFLG_JPLEPH:
      retc = jplplan(tjd, ipli, iflag, DO_SAVE, NULL, NULL, NULL, serr);

      if (retc == ERR)
	return ERR;

      if (retc == NOT_AVAILABLE) {
	iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_SWIEPH;
	if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	  strcat(serr, " \ntrying Swiss Eph; ");
	goto sweph_planet;
      } else if (retc == BEYOND_EPH_LIMITS) {
	if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
	  iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_MOSEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \nusing Moshier Eph; ");
	  goto moshier_planet;
	} else {
	  return ERR;
	}
      }

      if (ipli == SEI_SUN) {
	retc = app_pos_etc_sun(iflag, serr);
      } else {
	retc = app_pos_etc_plan(ipli, iplmoon, iflag, serr);
      }
      if (retc == ERR)
	return ERR;

      if (retc == NOT_AVAILABLE) {
	iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_SWIEPH;
	if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	  strcat(serr, " \ntrying Swiss Eph; ");
	goto sweph_planet;
      } else if (retc == BEYOND_EPH_LIMITS) {
	if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
	  iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_MOSEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \nusing Moshier Eph; ");
	  goto moshier_planet;
	} else
	  return ERR;
      }
      break;
    case SEFLG_SWIEPH:
      sweph_planet:

      retc = sweplan(tjd, ipli, SEI_FILE_PLANET, iflag, DO_SAVE, NULL, NULL, NULL, NULL, serr);
      if (retc == ERR)
	return ERR;

      if (retc == NOT_AVAILABLE) {
	if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
	  iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \nusing Moshier eph.; ");
	  goto moshier_planet;
	} else
	  return ERR;
      }

      if (ipli == SEI_SUN) {
	retc = app_pos_etc_sun(iflag, serr);
      } else {
	retc = app_pos_etc_plan(ipli, iplmoon, iflag, serr);
      }
      if (retc == ERR)
	return ERR;

      if (retc == NOT_AVAILABLE) {
	if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
	  iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \nusing Moshier eph.; ");
	  goto moshier_planet;
	} else
	  return ERR;
      }
      break;
    case SEFLG_MOSEPH:
      moshier_planet:
      retc = swi_moshplan(tjd, ipli, DO_SAVE, NULL, NULL, serr);
      if (retc == ERR)
	return ERR;

      if (ipli == SEI_SUN) {
	retc = app_pos_etc_sun(iflag, serr);
      } else {
	retc = app_pos_etc_plan(ipli, iplmoon, iflag, serr);
      }
      if (retc == ERR)
	return ERR;
      break;
    default:
      break;
  }
  return OK;
}

static int main_planet_bary(double tjd, int ipli, int32 epheflag, int32 iflag, AS_BOOL do_save,
		       double *xp, double *xe, double *xs, double *xm,
		       char *serr)
{
  int i, retc;
  switch(epheflag) {
    case SEFLG_JPLEPH:
      retc = jplplan(tjd, ipli, iflag, do_save, xp, xe, xs, serr);

      if (retc == ERR || retc == BEYOND_EPH_LIMITS)
	return retc;

      if (retc == NOT_AVAILABLE) {
	iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_SWIEPH;
	if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	  strcat(serr, " \ntrying Swiss Eph; ");
	goto sweph_planet;
      }
      break;
    case SEFLG_SWIEPH:
      sweph_planet:

      retc = sweplan(tjd, ipli, SEI_FILE_PLANET, iflag, do_save, xp, xe, xs, xm, serr);

      if (retc == ERR)
	return ERR;

      if (retc == NOT_AVAILABLE) {
	if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
	  iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \nusing Moshier eph.; ");
	  goto moshier_planet;
	} else {
	  return ERR;
	}
      }
      break;
    case SEFLG_MOSEPH:
      moshier_planet:
      retc = swi_moshplan(tjd, ipli, do_save, xp, xe, serr);
      if (retc == ERR)
	return ERR;
      for (i = 0; i <= 5; i++)
	xs[i] = 0;
      break;
    default:
      break;
  }
  return OK;
}

static int swemoon(double tjd, int32 iflag, AS_BOOL do_save, double *xpret, char *serr)
{
  int i, retc;
  struct plan_data *pdp = &swed.pldat[SEI_MOON];
  int32 speedf1, speedf2;
  double xx[6], *xp;
  if (do_save) {
    xp = pdp->x;
  } else {
    xp = xx;
  }

  speedf1 = pdp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == pdp->teval
	&& pdp->iephe == SEFLG_SWIEPH
	&& (!speedf2 || speedf1)) {
    xp = pdp->x;
  } else {

    retc = sweph(tjd, SEI_MOON, SEI_FILE_MOON, iflag, NULL, do_save, xp, serr);
    if (retc != OK)
      return(retc);
    if (do_save) {
      pdp->teval = tjd;
      pdp->xflgs = -1;
      pdp->iephe = SEFLG_SWIEPH;
    }
  }
  if (xpret != NULL)
    for (i = 0; i <= 5; i++)
      xpret[i] = xp[i];
  return(OK);
}

static int sweplan(double tjd, int ipli, int ifno, int32 iflag, AS_BOOL do_save,
		   double *xpret, double *xperet, double *xpsret, double *xpmret,
		   char *serr)
{
  int i, retc;
  int do_earth = FALSE, do_moon = FALSE, do_sunbary = FALSE;
  struct plan_data *pdp = &swed.pldat[ipli];
  struct plan_data *pebdp = &swed.pldat[SEI_EMB];
  struct plan_data *psbdp = &swed.pldat[SEI_SUNBARY];
  struct plan_data *pmdp = &swed.pldat[SEI_MOON];
  double xxp[6], xxm[6], xxs[6], xxe[6];
  double *xp, *xpe, *xpm, *xps;
  int32 speedf1, speedf2;

  if (do_save || ipli == SEI_SUNBARY || (pdp->iflg & SEI_FLG_HELIO)
    || xpsret != NULL || (iflag & SEFLG_HELCTR))
    do_sunbary = TRUE;
  if (do_save || ipli == SEI_EARTH || xperet != NULL)
    do_earth = TRUE;
  if (ipli == SEI_MOON) {
    do_earth = TRUE;
    do_sunbary = TRUE;
  }
  if (do_save || ipli == SEI_MOON || ipli == SEI_EARTH || xperet != NULL || xpmret != NULL)
    do_moon = TRUE;
  if (do_save)  {
    xp = pdp->x;
    xpe = pebdp->x;
    xps = psbdp->x;
    xpm = pmdp->x;
  } else {
    xp = xxp;
    xpe = xxe;
    xps = xxs;
    xpm = xxm;
  }
  speedf2 = iflag & SEFLG_SPEED;

  if (do_sunbary) {
    speedf1 = psbdp->xflgs & SEFLG_SPEED;

    if (tjd == psbdp->teval
	  && psbdp->iephe == SEFLG_SWIEPH
	  && (!speedf2 || speedf1)) {
      for (i = 0; i <= 5; i++)
	xps[i] = psbdp->x[i];
    } else {
      retc = sweph(tjd, SEI_SUNBARY, SEI_FILE_PLANET, iflag, NULL, do_save, xps, serr);
      if (retc != OK)
	return(retc);
    }
    if (xpsret != NULL)
      for (i = 0; i <= 5; i++)
	xpsret[i] = xps[i];
  }

  if (do_moon) {
    speedf1 = pmdp->xflgs & SEFLG_SPEED;
    if (tjd == pmdp->teval
	  && pmdp->iephe == SEFLG_SWIEPH
	  && (!speedf2 || speedf1)) {
      for (i = 0; i <= 5; i++)
	xpm[i] = pmdp->x[i];
    } else {
      retc = sweph(tjd, SEI_MOON, SEI_FILE_MOON, iflag, NULL, do_save, xpm, serr);
      if (retc == ERR)
	return(retc);

      if (swed.fidat[SEI_FILE_MOON].fptr == NULL) {
	if (serr != NULL && strlen(serr) + 35 < AS_MAXCH)
	  strcat(serr, " \nusing Moshier eph. for moon; ");
	retc = swi_moshmoon(tjd, do_save, xpm, serr);
	if (retc != OK)
	  return(retc);
      }
    }
    if (xpmret != NULL)
      for (i = 0; i <= 5; i++)
	xpmret[i] = xpm[i];
  }

  if (do_earth) {
    speedf1 = pebdp->xflgs & SEFLG_SPEED;
    if (tjd == pebdp->teval
	  && pebdp->iephe == SEFLG_SWIEPH
	  && (!speedf2 || speedf1)) {
      for (i = 0; i <= 5; i++)
	xpe[i] = pebdp->x[i];
    } else {
      retc = sweph(tjd, SEI_EMB, SEI_FILE_PLANET, iflag, NULL, do_save, xpe, serr);
      if (retc != OK)
	return(retc);

      embofs(xpe, xpm);

      if (xpe == pebdp->x || (iflag & SEFLG_SPEED))
	embofs(xpe+3, xpm+3);
    }
    if (xperet != NULL)
      for (i = 0; i <= 5; i++)
	xperet[i] = xpe[i];
  }
  if (ipli == SEI_MOON) {
    for (i = 0; i <= 5; i++)
      xp[i] = xpm[i];
  } else if (ipli == SEI_EARTH) {
    for (i = 0; i <= 5; i++)
      xp[i] = xpe[i];
  } else if (ipli == SEI_SUN) {
    for (i = 0; i <= 5; i++)
      xp[i] = xps[i];
  } else {

    speedf1 = pdp->xflgs & SEFLG_SPEED;
    if (tjd == pdp->teval
	  && pdp->iephe == SEFLG_SWIEPH
	  && (!speedf2 || speedf1)) {
      for (i = 0; i <= 5; i++)
	xp[i] = pdp->x[i];
      return(OK);
    } else {
      retc = sweph(tjd, ipli, ifno, iflag, NULL, do_save, xp, serr);
      if (retc != OK)
	return(retc);

      if (pdp->iflg & SEI_FLG_HELIO) {

	for (i = 0; i <= 2; i++)
	  xp[i] += xps[i];
	if (do_save || (iflag & SEFLG_SPEED))
	  for (i = 3; i <= 5; i++)
	    xp[i] += xps[i];
      }
    }
  }
  if (xpret != NULL)
    for (i = 0; i <= 5; i++)
      xpret[i] = xp[i];
  return(OK);
}

static int jplplan(double tjd, int ipli, int32 iflag, AS_BOOL do_save,
		   double *xpret, double *xperet, double *xpsret, char *serr)
{
  int i, retc;
  AS_BOOL do_earth = FALSE, do_sunbary = FALSE;
  double ss[3];
  double xxp[6], xxe[6], xxs[6];
  double *xp, *xpe, *xps;
  int ictr = J_SBARY;
  struct plan_data *pdp = &swed.pldat[ipli];
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  iflag = SEFLG_JPLEPH;

  if (do_save) {
    xp = pdp->x;
    xpe = pedp->x;
    xps = psdp->x;
  } else {
    xp = xxp;
    xpe = xxe;
    xps = xxs;
  }
  if (do_save || ipli == SEI_EARTH || xperet != NULL
    || (ipli == SEI_MOON))
    do_earth = TRUE;
  if (do_save || ipli == SEI_SUNBARY || xpsret != NULL
    || (ipli == SEI_MOON))
    do_sunbary = TRUE;
  if (ipli == SEI_MOON)
    ictr = J_EARTH;

  if (!swed.jpl_file_is_open) {
    retc = open_jpl_file(ss, swed.jplfnam, swed.ephepath, serr);
    if (retc != OK)
      return (retc);
  }
  if (do_earth) {

    if (tjd != pedp->teval || tjd == 0) {
      retc = swi_pleph(tjd, J_EARTH, J_SBARY, xpe, serr);
      if (do_save) {
	pedp->teval = tjd;
	pedp->xflgs = -1;
	pedp->iephe = SEFLG_JPLEPH;
      }
      if (retc != OK) {
	swi_close_jpl_file();
	swed.jpl_file_is_open = FALSE;
	return retc;
      }
    } else {
      xpe = pedp->x;
    }
    if (xperet != NULL)
      for (i = 0; i <= 5; i++)
	xperet[i] = xpe[i];

  }
  if (do_sunbary) {

    if (tjd != psdp->teval || tjd == 0) {
      retc = swi_pleph(tjd, J_SUN, J_SBARY, xps, serr);
      if (do_save) {
	psdp->teval = tjd;
	psdp->xflgs = -1;
	psdp->iephe = SEFLG_JPLEPH;
      }
      if (retc != OK) {
	swi_close_jpl_file();
	swed.jpl_file_is_open = FALSE;
	return retc;
      }
    } else {
      xps = psdp->x;
    }
    if (xpsret != NULL)
      for (i = 0; i <= 5; i++)
	xpsret[i] = xps[i];
  }

  if (ipli == SEI_EARTH) {
    for (i = 0; i <= 5; i++)
      xp[i] = xpe[i];

  } if (ipli == SEI_SUNBARY) {
    for (i = 0; i <= 5; i++)
      xp[i] = xps[i];

  } else {

    if (tjd == pdp->teval && pdp->iephe == SEFLG_JPLEPH) {
      xp = pdp->x;
    } else {
      retc = swi_pleph(tjd, pnoint2jpl[ipli], ictr, xp, serr);
      if (do_save) {
	pdp->teval = tjd;
	pdp->xflgs = -1;
	pdp->iephe = SEFLG_JPLEPH;
      }
      if (retc != OK) {
	swi_close_jpl_file();
	swed.jpl_file_is_open = FALSE;
	return retc;
      }
    }
  }
  if (xpret != NULL)
    for (i = 0; i <= 5; i++)
      xpret[i] = xp[i];
  return (OK);
}

static int sweph(double tjd, int ipli, int ifno, int32 iflag, double *xsunb, AS_BOOL do_save, double *xpret, char *serr)
{
  int i, ipl, retc, subdirlen;
  char s[2 * AS_MAXCH], subdirnam[AS_MAXCH], fname[AS_MAXCH], *sp;
  double t, tsv;
  double xemb[6], xx[6], *xp;
  struct plan_data *pdp;
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  struct file_data *fdp = &swed.fidat[ifno];
  int32 speedf1, speedf2;
  AS_BOOL need_speed;
  ipl = ipli;
  if (ipli > SE_AST_OFFSET)
    ipl = SEI_ANYBODY;
  if (ipli > SE_PLMOON_OFFSET)
    ipl = SEI_ANYBODY;
  pdp = &swed.pldat[ipl];
  if (do_save) {
    xp = pdp->x;
  } else {
    xp = xx;
  }

  speedf1 = pdp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == pdp->teval
	&& pdp->iephe == SEFLG_SWIEPH
	&& (!speedf2 || speedf1)
        && ipl < SEI_ANYBODY) {
    if (xpret != NULL)
      for (i = 0; i <= 5; i++)
	xpret[i] = pdp->x[i];
    return(OK);
  }

  if (fdp->fptr != NULL) {

    if (tjd < fdp->tfstart || tjd > fdp->tfend
      || (ipl == SEI_ANYBODY && ipli != pdp->ibdy)) {
      fclose(fdp->fptr);
      fdp->fptr = NULL;
      if (pdp->refep != NULL)
	free((void *) pdp->refep);
      pdp->refep = NULL;
      if (pdp->segp != NULL)
	free((void *) pdp->segp);
      pdp->segp = NULL;
    }
  }

  if (fdp->fptr == NULL) {
    swi_gen_filename(tjd, ipli, fname);
    strcpy(subdirnam, fname);
    sp = strrchr(subdirnam, (int) *DIR_GLUE);
    if (sp != NULL) {
      *sp = '\0';
      subdirlen = (int) strlen(subdirnam);
    } else {
      subdirlen = 0;
    }
    strcpy(s, fname);
again:
    fdp->fptr = swi_fopen(ifno, s, swed.ephepath, serr);
    if (fdp->fptr == NULL) {

      if (ipli > SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET) {
	if (subdirlen > 0 && strncmp(s, subdirnam, (size_t) subdirlen) == 0) {
	  swi_strcpy(s, s + subdirlen + 1);
	  goto again;
	}

      } else if (ipli > SE_AST_OFFSET) {
	char *spp;
	spp = strchr(s, '.');
	if (spp > s && *(spp-1) != 's') {
	  sprintf(spp, "s.%s", SE_FILE_SUFFIX);
	  goto again;
	}

        spp--;
	swi_strcpy(spp, spp + 1);
	if (subdirlen > 0 && strncmp(s, subdirnam, (size_t) subdirlen) == 0) {
	  swi_strcpy(s, s + subdirlen + 1);
	  goto again;
	}
      }
      return(NOT_AVAILABLE);
    }

    if (serr != NULL) *serr = '\0';
    retc = read_const(ifno, serr);
    if (retc != OK)
      return(retc);
  }

  if (tjd < fdp->tfstart || tjd > fdp->tfend) {
    if (serr != NULL) {
      sp = strrchr(fname, (int) *DIR_GLUE);
      if (sp != NULL) {
        sp++;
      } else {
        sp = fname;
      }
      if (ipli > SE_AST_OFFSET) {
        sprintf(s, "asteroid No. %d (%s): ", ipli - SE_AST_OFFSET, sp);
      } else if (ipli > SE_PLMOON_OFFSET) {
	if (strstr(fname, "99.") != NULL)
	  sprintf(s, "plan. COB No. %d (%s): ", ipli, sp);
	else
	  sprintf(s, "plan. moon No. %d (%s): ", ipli, sp);
      } else if (ipli > SEI_PLUTO) {
        sprintf(s, "asteroid eph. file (%s): ", sp);
      } else if (ipli != SEI_MOON) {
        sprintf(s, "planets eph. file (%s): ", sp);
      } else {
        sprintf(s, "moon eph. file (%s): ", sp);
      }
      if (tjd < fdp->tfstart) {
	sprintf(s + strlen(s), "jd %f < lower limit %f;",
		  tjd, fdp->tfstart);
      } else {
	sprintf(s + strlen(s), "jd %f > upper limit %f;",
		  tjd, fdp->tfend);
      }
      if (strlen(serr) + strlen(s) < AS_MAXCH)
	strcat(serr, s);
    }
    return(NOT_AVAILABLE);
  }

  if (pdp->segp == NULL || tjd < pdp->tseg0 || tjd > pdp->tseg1) {
    retc = get_new_segment(tjd, ipl, ifno, serr);
    if (retc != OK)
      return(retc);

    if (pdp->iflg & SEI_FLG_ROTATE) {
      rot_back(ipl);
    } else {
      pdp->neval = pdp->ncoe;
    }
  }

  t = (tjd - pdp->tseg0) / pdp->dseg;
  t = t * 2 - 1;

  need_speed = (do_save || (iflag & SEFLG_SPEED));
  for (i = 0; i <= 2; i++) {
    xp[i]  = swi_echeb (t, pdp->segp+(i*pdp->ncoe), pdp->neval);
    if (need_speed) {
      xp[i+3] = swi_edcheb(t, pdp->segp+(i*pdp->ncoe), pdp->neval) / pdp->dseg * 2;
    } else {
      xp[i+3] = 0;
    }
  }

  if (ipl == SEI_SUNBARY && (pdp->iflg & SEI_FLG_EMBHEL)) {

    tsv = pedp->teval;
    pedp->teval = 0;
    retc = sweph(tjd, SEI_EMB, ifno, iflag | SEFLG_SPEED, NULL, NO_SAVE, xemb, serr);
    if (retc != OK)
      return(retc);
    pedp->teval = tsv;
    for (i = 0; i <= 2; i++)
      xp[i] = xemb[i] - xp[i];
    if (need_speed)
      for (i = 3; i <= 5; i++)
	xp[i] = xemb[i] - xp[i];
  }
#if 1

  if (xsunb != NULL && ((iflag & SEFLG_JPLEPH) || (iflag & SEFLG_SWIEPH))) {
    if (ipl >= SEI_ANYBODY) {
      for (i = 0; i <= 2; i++)
	xp[i] += xsunb[i];
      if (need_speed)
	for (i = 3; i <= 5; i++)
	  xp[i] += xsunb[i];
    }
  }
#endif
  if (do_save) {
    pdp->teval = tjd;
    pdp->xflgs = -1;
    if (ifno == SEI_FILE_PLANET || ifno == SEI_FILE_MOON) {
      pdp->iephe = SEFLG_SWIEPH;
    } else {
      pdp->iephe = psdp->iephe;
    }
  }
  if (xpret != NULL)
    for (i = 0; i <= 5; i++)
      xpret[i] = xp[i];
  return(OK);
}

FILE *swi_fopen(int ifno, char *fname, char *ephepath, char *serr)
{
  int np, i, j;
  FILE *fp = NULL;
  char *fnamp, fn[AS_MAXCH];
  char *cpos[20];
  char s[2 * AS_MAXCH];
  char s1[AS_MAXCH];
  if (ifno >= 0) {
    fnamp = swed.fidat[ifno].fnam;
  } else {
    fnamp = fn;
  }
  strcpy(s1, ephepath);
  np = swi_cutstr(s1, PATH_SEPARATOR, cpos, 20);
  *s = '\0';
  for (i = 0; i < np; i++) {
    strcpy(s, cpos[i]);
    if (strcmp(s, ".") == 0) {
      *s = '\0';
    } else {
      j = (int) strlen(s);
      if (*s != '\0' && *(s + j - 1) != *DIR_GLUE)
	strcat(s, DIR_GLUE);
    }
    if (strlen(s) + strlen(fname) < AS_MAXCH) {
      strcat(s, fname);
    } else {
      if (serr != NULL)
	sprintf(serr, "error: file path and name must be shorter than %d.", AS_MAXCH);
      return NULL;
    }
    strcpy(fnamp, s);
    fp = fopen(fnamp, BFILE_R_ACCESS);
    if (fp != NULL)
      return fp;
  }
  sprintf(s, "SwissEph file '%s' not found in PATH '%s'", fname, ephepath);
  s[AS_MAXCH-1] = '\0';
  if (serr != NULL)
    strcpy(serr, s);
  return NULL;
}

int32 swi_get_denum(int32 ipli, int32 iflag)
{
  struct file_data *fdp = NULL;
  if (iflag & SEFLG_MOSEPH)
    return 403;
  if (iflag & SEFLG_JPLEPH) {
    if (swed.jpldenum > 0) {
      return swed.jpldenum;
    } else {
      return SE_DE_NUMBER;
    }
  }
  if (ipli > SE_AST_OFFSET) {
    fdp = &swed.fidat[SEI_FILE_ANY_AST];
  } else if (ipli > SE_PLMOON_OFFSET) {
    fdp = &swed.fidat[SEI_FILE_ANY_AST];
  } else if (ipli == SEI_CHIRON
      || ipli == SEI_PHOLUS
      || ipli == SEI_CERES
      || ipli == SEI_PALLAS
      || ipli == SEI_JUNO
      || ipli == SEI_VESTA) {
    fdp = &swed.fidat[SEI_FILE_MAIN_AST];
  } else if (ipli == SEI_MOON) {
    fdp = &swed.fidat[SEI_FILE_MOON];
  } else {
    fdp = &swed.fidat[SEI_FILE_PLANET];
  }
  if (fdp != NULL) {
    if (fdp->sweph_denum != 0) {
      return fdp->sweph_denum;
    } else {
      return SE_DE_NUMBER;
    }
  }
  return SE_DE_NUMBER;
}

static int calc_center_body(int32 ipli, int32 iflag, double *xx, double *xcom, char *serr)
{
  int i;
  if (!(iflag & SEFLG_CENTER_BODY))
    return OK;
  if (ipli < SEI_MARS || ipli > SEI_PLUTO)
    return OK;
  for (i = 0; i <= 5; i++)
    xx[i] += xcom[i];
  return OK;
}

static int app_pos_etc_plan(int ipli, int iplmoon, int32 iflag, char *serr)
{
  int i, j, niter, retc = OK;
  int ipl, ifno, ibody;
  int32 flg1, flg2;
  double xx[6], xx0[6], dx[3], dt, t, dtsave_for_defl;
  double xobs[6], xobs2[6];
  double xearth[6], xsun[6], xcom[6];
  double xxsp[6], xxsv[6];
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *pdp;
  struct epsilon *oe = &swed.oec2000;
  int32 epheflag = iflag & SEFLG_EPHMASK;
  dtsave_for_defl = 0;

  if (ipli > SE_PLMOON_OFFSET || ipli > SE_AST_OFFSET) {
    ifno = SEI_FILE_ANY_AST;
    ibody = IS_ANY_BODY;
    pdp = &swed.pldat[SEI_ANYBODY];
  } else if (ipli == SEI_CHIRON
      || ipli == SEI_PHOLUS
      || ipli == SEI_CERES
      || ipli == SEI_PALLAS
      || ipli == SEI_JUNO
      || ipli == SEI_VESTA) {
    ifno = SEI_FILE_MAIN_AST;
    ibody = IS_MAIN_ASTEROID;
    pdp = &swed.pldat[ipli];
  } else {
    ifno = SEI_FILE_PLANET;
    ibody = IS_PLANET;
    pdp = &swed.pldat[ipli];
  }
  t = pdp->teval;

  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pdp->xflgs = iflag;
    pdp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }

  for (i = 0; i <= 5; i++)
    xx[i] = pdp->x[i];

  calc_center_body(ipli, iflag, xx, swed.pldat[SEI_ANYBODY].x, serr);
  for (i = 0; i <= 5; i++)
    xx0[i] = xx[i];

  if (iflag & SEFLG_HELCTR) {
    if (pdp->iephe == SEFLG_JPLEPH || pdp->iephe == SEFLG_SWIEPH)
      for (i = 0; i <= 5; i++)
	xx[i] -= swed.pldat[SEI_SUNBARY].x[i];
  }

  if (iflag & SEFLG_TOPOCTR) {
    if (swed.topd.teval != pedp->teval
      || swed.topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr) != OK)
        return ERR;
    } else {
      for (i = 0; i <= 5; i++)
        xobs[i] = swed.topd.xobs[i];
    }

    for (i = 0; i <= 5; i++)
      xobs[i] = xobs[i] + pedp->x[i];
  } else {

    for (i = 0; i <= 5; i++)
      xobs[i] = pedp->x[i];
  }

  if (!(iflag & SEFLG_TRUEPOS)) {

    if (pdp->iephe == SEFLG_JPLEPH || pdp->iephe == SEFLG_SWIEPH) {
      niter = 1;
    } else {
      niter = 0;
    }
    if (iflag & SEFLG_SPEED) {

      for (i = 0; i <= 2; i++)
	xxsv[i] = xxsp[i] = xx[i] - xx[i+3];
      for (j = 0; j <= niter; j++) {
	for (i = 0; i <= 2; i++) {
	  dx[i] = xxsp[i];
	  if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR))
	    dx[i] -= (xobs[i] - xobs[i+3]);
	}

	dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
	for (i = 0; i <= 2; i++) {

	  xxsp[i] = xxsv[i] - dt * xx0[i+3];
	}
      }

      for (i = 0; i <= 2; i++)
	xxsp[i] = xxsv[i] - xxsp[i];
    }

    for (j = 0; j <= niter; j++) {
      for (i = 0; i <= 2; i++) {
	dx[i] = xx[i];
	if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR))
	  dx[i] -= xobs[i];
      }
      dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;

      t = pdp->teval - dt;
      dtsave_for_defl = dt;
      for (i = 0; i <= 2; i++) {

	xx[i] = xx0[i] - dt * xx0[i+3];
      }
    }

    if (iflag & SEFLG_SPEED) {
      for (i = 0; i <= 2; i++) {

	xxsp[i] = xx0[i] - xx[i] - xxsp[i];
      }
    }

    if ((iflag & SEFLG_CENTER_BODY)
      && ipli >= SE_MARS && ipli <= SE_PLUTO) {

      retc = sweph(t, iplmoon, SEI_FILE_ANY_AST, iflag, NULL, NO_SAVE, xcom, serr);
      if (retc == ERR || retc == NOT_AVAILABLE)
	return ERR;
    }
    switch(epheflag) {
      case SEFLG_JPLEPH:
	if (ibody >= IS_ANY_BODY)
	  ipl = -1;
	else
	  ipl = pnoint2jpl[ipli];
	if (ibody == IS_PLANET) {
	  retc = swi_pleph(t, ipl, J_SBARY, xx, serr);
	  if (retc != OK) {
	    swi_close_jpl_file();
	    swed.jpl_file_is_open = FALSE;
	  }
	} else {

	  retc = swi_pleph(t, J_SUN, J_SBARY, xsun, serr);
	  if (retc != OK) {
	    swi_close_jpl_file();
	    swed.jpl_file_is_open = FALSE;
	  }

	  retc = sweph(t, ipli, ifno, iflag, xsun, NO_SAVE, xx, serr);
	}
	if (retc != OK)
	  return(retc);

	if ((iflag & SEFLG_SPEED)
	  && !(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR)) {
	  retc = swi_pleph(t, J_EARTH, J_SBARY, xearth, serr);
	  if (retc != OK) {
	    swi_close_jpl_file();
	    swed.jpl_file_is_open = FALSE;
	    return(retc);
	  }
	}
	break;
      case SEFLG_SWIEPH:
	if (ibody == IS_PLANET) {
	  retc = sweplan(t, ipli, ifno, iflag, NO_SAVE, xx, xearth, xsun, NULL, serr);
	} else {
	  retc = sweplan(t, SEI_EARTH, SEI_FILE_PLANET, iflag, NO_SAVE, xearth, NULL, xsun, NULL, serr);
	  if (retc == OK)
	    retc = sweph(t, ipli, ifno, iflag, xsun, NO_SAVE, xx, serr);
	}
	if (retc != OK)
	  return(retc);
	break;
      case SEFLG_MOSEPH:
      default:

	if (iflag & SEFLG_SPEED
	  && !(iflag & (SEFLG_HELCTR | SEFLG_BARYCTR))) {
	  if (ibody == IS_PLANET) {
	    retc = swi_moshplan(t, ipli, NO_SAVE, xxsv, xearth, serr);
          } else {
	    retc = sweph(t, ipli, ifno, iflag, NULL, NO_SAVE, xxsv, serr);
	    if (retc == OK)
	      retc = swi_moshplan(t, SEI_EARTH, NO_SAVE, xearth, xearth, serr);
          }
	  if (retc != OK)
	    return(retc);

	  for (i = 3; i <= 5; i++)
	    xx[i] = xxsv[i];
        }
	break;
    }
    calc_center_body(ipli, iflag, xx, xcom, serr);
    if (iflag & SEFLG_HELCTR) {
      if (pdp->iephe == SEFLG_JPLEPH || pdp->iephe == SEFLG_SWIEPH)
	for (i = 0; i <= 5; i++)
	  xx[i] -= swed.pldat[SEI_SUNBARY].x[i];
    }
    if (iflag & SEFLG_SPEED) {

      if (iflag & SEFLG_TOPOCTR) {
        if (swi_get_observer(t, iflag | SEFLG_NONUT, NO_SAVE, xobs2, serr) != OK)
          return ERR;
        for (i = 0; i <= 5; i++)
          xobs2[i] += xearth[i];
      } else {
        for (i = 0; i <= 5; i++)
          xobs2[i] = xearth[i];
      }
    }
  }

  if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR)) {

    for (i = 0; i <= 5; i++)
      xx[i] -= xobs[i];
    if ((iflag & SEFLG_TRUEPOS) == 0 ) {

      if (iflag & SEFLG_SPEED)
	for (i = 3; i <= 5; i++)
	  xx[i] -= xxsp[i-3];
    }
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOGDEFL))

    swi_deflect_light(xx, dtsave_for_defl, iflag);

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOABERR)) {

    swi_aberr_light(xx, xobs, iflag);

    if (iflag & SEFLG_SPEED) {
      for (i = 3; i <= 5; i++)
	xx[i] += xobs[i] - xobs2[i];
    }
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_ICRS) && swi_get_denum(ipli, epheflag) >= 403) {
    swi_bias(xx, t, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];

  if (!(iflag & SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, pdp->teval, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else {
    oe = &swed.oec2000;
  }
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr);
}

static int app_pos_rest(struct plan_data *pdp, int32 iflag,
                        double *xx, double *x2000,
                        struct epsilon *oe, char *serr)
{
  int i;
  double daya[2];
  double xxsv[24];

  if (!(iflag & SEFLG_NONUT))
    swi_nutate(xx, iflag, FALSE);

  for (i = 0; i <= 5; i++)
    pdp->xreturn[18+i] = xx[i];

  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  if (iflag & SEFLG_SPEED)
    swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
  if (!(iflag & SEFLG_NONUT)) {
    swi_coortrf2(xx, xx, swed.nut.snut, swed.nut.cnut);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(xx+3, xx+3, swed.nut.snut, swed.nut.cnut);
  }

  for (i = 0; i <= 5; i++)
    pdp->xreturn[6+i] = xx[i];

  if (iflag & SEFLG_SIDEREAL) {

    if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
      if (swi_trop_ra2sid_lon(x2000, pdp->xreturn+6, pdp->xreturn+18, iflag) != OK)
	return ERR;

    } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
      if (swi_trop_ra2sid_lon_sosy(x2000, pdp->xreturn+6, iflag) != OK)
	return ERR;
    } else {

      swi_cartpol_sp(pdp->xreturn+6, pdp->xreturn);

      for (i = 0; i < 24; i++)
        xxsv[i] = pdp->xreturn[i];
      if (swi_get_ayanamsa_with_speed(pdp->teval, iflag, daya, serr) == ERR)
        return ERR;

      for (i = 0; i < 24; i++)
        pdp->xreturn[i] = xxsv[i];
      pdp->xreturn[0] -= daya[0] * DEGTORAD;
      pdp->xreturn[3] -= daya[1] * DEGTORAD;
      swi_polcart_sp(pdp->xreturn, pdp->xreturn+6);
    }
  }

  swi_cartpol_sp(pdp->xreturn+18, pdp->xreturn+12);
  swi_cartpol_sp(pdp->xreturn+6, pdp->xreturn);

    for (i = 0; i < 2; i++) {
      pdp->xreturn[i] *= RADTODEG;
      pdp->xreturn[i+3] *= RADTODEG;
      pdp->xreturn[i+12] *= RADTODEG;
      pdp->xreturn[i+15] *= RADTODEG;
    }

  pdp->xflgs = iflag;
  pdp->iephe = iflag & SEFLG_EPHMASK;
  return OK;
}

void CALL_CONV swe_set_sid_mode(int32 sid_mode, double t0, double ayan_t0)
{
  struct sid_data *sip = &swed.sidd;
  swi_init_swed_if_start();
  if (sid_mode < 0)
    sid_mode = 0;
  sip->sid_mode = sid_mode;
  if (sid_mode >= SE_SIDBITS)
    sid_mode %= SE_SIDBITS;

  if (sid_mode == SE_SIDM_J2000
	  || sid_mode == SE_SIDM_J1900
	  || sid_mode == SE_SIDM_B1950
	  || sid_mode == SE_SIDM_GALALIGN_MARDYKS
	  ) {

    sip->sid_mode = sid_mode;
    sip->sid_mode |= SE_SIDBIT_ECL_T0;
  }
  if (sid_mode == SE_SIDM_TRUE_CITRA
      || sid_mode == SE_SIDM_TRUE_REVATI
      || sid_mode == SE_SIDM_TRUE_PUSHYA
      || sid_mode == SE_SIDM_TRUE_SHEORAN
      || sid_mode == SE_SIDM_TRUE_MULA
      || sid_mode == SE_SIDM_GALCENT_0SAG
      || sid_mode == SE_SIDM_GALCENT_COCHRANE
      || sid_mode == SE_SIDM_GALCENT_RGILBRAND
      || sid_mode == SE_SIDM_GALCENT_MULA_WILHELM
      || sid_mode == SE_SIDM_GALEQU_IAU1958
      || sid_mode == SE_SIDM_GALEQU_TRUE
      || sid_mode == SE_SIDM_GALEQU_MULA
      ) {

    sip->sid_mode = sid_mode;
  }

  if (sid_mode >= SE_NSIDM_PREDEF && sid_mode != SE_SIDM_USER)
    sip->sid_mode = sid_mode = SE_SIDM_FAGAN_BRADLEY;
  swed.ayana_is_set = TRUE;
  if (sid_mode == SE_SIDM_USER) {
    sip->t0 = t0;
    sip->ayan_t0 = ayan_t0;
    sip->t0_is_UT = FALSE;
    if (sip->sid_mode & SE_SIDBIT_USER_UT)
      sip->t0_is_UT = TRUE;
  } else {
    sip->t0 = ayanamsa[sid_mode].t0;
    sip->ayan_t0 = ayanamsa[sid_mode].ayan_t0;
    sip->t0_is_UT = ayanamsa[sid_mode].t0_is_UT;
  }

  if (sid_mode < SE_NSIDM_PREDEF && (sip->sid_mode & SE_SIDBIT_PREC_ORIG) && ayanamsa[sid_mode].prec_offset > 0) {
    swed.astro_models[SE_MODEL_PREC_LONGTERM] = ayanamsa[sid_mode].prec_offset;
    swed.astro_models[SE_MODEL_PREC_SHORTTERM] = ayanamsa[sid_mode].prec_offset;

    switch(ayanamsa[sid_mode].prec_offset) {
      case SEMOD_PREC_NEWCOMB:
        swed.astro_models[SE_MODEL_NUT] = SEMOD_NUT_WOOLARD;
	break;
      case SEMOD_PREC_IAU_1976:
        swed.astro_models[SE_MODEL_NUT] = SEMOD_NUT_IAU_1980;
	break;
      default:
        break;
    }
  }
  swi_force_app_pos_etc();
}

int32 CALL_CONV swe_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr)
{
  struct nut nuttmp;
  struct nut *nutp = &nuttmp;
  int32 retval = swi_get_ayanamsa_ex(tjd_et, iflag, daya, serr);
  if (!(iflag & SEFLG_NONUT)) {
    if (tjd_et == swed.nut.tnut) {
      nutp = &swed.nut;
    } else {
      nutp = &nuttmp;
      swi_nutation(tjd_et, iflag, nutp->nutlo);
    }
    *daya += nutp->nutlo[0] * RADTODEG;
    retval &= (~SEFLG_NONUT);
  }
  return retval;
}

static int get_aya_correction(int iflag, double *corr, char *serr) {
  double x[6], eps, t0;
  struct sid_data *sip = &swed.sidd;
  int prec_model = swed.astro_models[SE_MODEL_PREC_LONGTERM];
  int prec_model_short = swed.astro_models[SE_MODEL_PREC_SHORTTERM];
  int prec_offset = 0;
  int sid_mode = sip->sid_mode;
  sid_mode %= SE_SIDBITS;
  *corr = 0;
  if (sip->t0 == J2000)
    return 0;
  if (sip->sid_mode & SE_SIDBIT_NO_PREC_OFFSET)
    return 0;
  if (sid_mode < SE_NSIDM_PREDEF)
    prec_offset = ayanamsa[sid_mode].prec_offset;
  if (prec_offset < 0) prec_offset = 0;
  if (prec_model == prec_offset)
    return 0;
  t0 = sip->t0;
  if (sip->t0_is_UT)
    t0 += swe_deltat_ex(t0, iflag, serr);

  x[0] = 1;
  x[1] = x[2] = 0;
  swi_precess(x, t0, 0, J_TO_J2000);
  swed.astro_models[SE_MODEL_PREC_LONGTERM] = prec_offset;
  swed.astro_models[SE_MODEL_PREC_SHORTTERM] = prec_offset;
  swi_precess(x, t0, 0, J2000_TO_J);
  swed.astro_models[SE_MODEL_PREC_LONGTERM] = prec_model;
  swed.astro_models[SE_MODEL_PREC_SHORTTERM] = prec_model_short;

  eps = swi_epsiln(t0, 0);
  swi_coortrf(x, x, eps);

  swi_cartpol(x, x);

  *corr = x[0] * RADTODEG;
  if (*corr > 350 ) *corr -= 360;

  return OK;
}

int32 swi_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr)
{
  double x[6], eps, t0, corr;
  struct sid_data *sip = &swed.sidd;
  char star[AS_MAXCH];
  int32 epheflag, otherflag, retflag, iflag_true, iflag_galequ;
  int sid_mode = sip->sid_mode;
  iflag = plaus_iflag(iflag, -1, tjd_et, serr);
  epheflag = iflag & SEFLG_EPHMASK;
  otherflag = iflag & ~SEFLG_EPHMASK;
  *daya = 0.0;
  iflag &= SEFLG_EPHMASK;
  iflag |= SEFLG_NONUT;
  sid_mode %= SE_SIDBITS;

  iflag_galequ = iflag | SEFLG_TRUEPOS;
#if 1

  iflag_true = iflag;
  if (otherflag & SEFLG_TRUEPOS) iflag_true |= SEFLG_TRUEPOS;
  if (otherflag & SEFLG_NOABERR) iflag_true |= SEFLG_NOABERR;
  if (otherflag & SEFLG_NOGDEFL) iflag_true |= SEFLG_NOGDEFL;
#endif

  if (swi_init_swed_if_start() == 1 && !(epheflag & SEFLG_MOSEPH)
     && (sid_mode ==  SE_SIDM_TRUE_CITRA
      || sid_mode == SE_SIDM_TRUE_REVATI
      || sid_mode == SE_SIDM_TRUE_PUSHYA
      || sip->sid_mode == SE_SIDM_TRUE_SHEORAN
      || sid_mode == SE_SIDM_TRUE_MULA
      || sid_mode == SE_SIDM_GALCENT_0SAG
      || sid_mode == SE_SIDM_GALCENT_COCHRANE
      || sid_mode == SE_SIDM_GALCENT_RGILBRAND
      || sid_mode == SE_SIDM_GALCENT_MULA_WILHELM
      || sid_mode == SE_SIDM_GALEQU_IAU1958
      || sid_mode == SE_SIDM_GALEQU_TRUE
      || sid_mode == SE_SIDM_GALEQU_MULA)
      && serr != NULL) {
    strcpy(serr, "Please call swe_set_ephe_path() or swe_set_jplfile() before calling swe_get_ayanamsa_ex()");
  }
  if (!swed.ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);
  if (sid_mode == SE_SIDM_TRUE_CITRA) {
    strcpy(star, "Spica");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR) {
      return ERR;
    }

    *daya = swe_degnorm(x[0] - 180);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode == SE_SIDM_TRUE_REVATI) {
    strcpy(star, ",zePsc");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 359.8333333333);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode == SE_SIDM_TRUE_PUSHYA) {
    strcpy(star, ",deCnc");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 106);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode == SE_SIDM_TRUE_SHEORAN) {
    strcpy(star, ",deCnc");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 103.49264221625);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode == SE_SIDM_TRUE_MULA) {
    strcpy(star, ",laSco");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 240);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode ==  SE_SIDM_GALCENT_0SAG) {
    strcpy(star, ",SgrA*");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 240.0);
    return (retflag & SEFLG_EPHMASK);

  }
  if (sid_mode ==  SE_SIDM_GALCENT_COCHRANE) {
    strcpy(star, ",SgrA*");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 270.0);
    return (retflag & SEFLG_EPHMASK);

  }
  if (sid_mode ==  SE_SIDM_GALCENT_RGILBRAND) {
    strcpy(star, ",SgrA*");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_true, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 210.0 - 90.0 * 0.3819660113);
    return (retflag & SEFLG_EPHMASK);

  }
  if (sid_mode == SE_SIDM_GALCENT_MULA_WILHELM) {
    strcpy(star, ",SgrA*");

    if ((retflag = swe_fixstar(star, tjd_et, iflag_true | SEFLG_EQUATORIAL, x, serr)) == ERR)
      return ERR;
    eps = swi_epsiln(tjd_et, iflag) * RADTODEG;
    *daya = swi_armc_to_mc(x[0], eps);
    *daya = swe_degnorm(*daya - 246.6666666667);
    return (retflag & SEFLG_EPHMASK);

  }
  if (sid_mode == SE_SIDM_GALEQU_IAU1958) {
    strcpy(star, ",GP1958");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_galequ, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 150);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode == SE_SIDM_GALEQU_TRUE) {
    strcpy(star, ",GPol");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_galequ, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 150);
    return (retflag & SEFLG_EPHMASK);
  }
  if (sid_mode == SE_SIDM_GALEQU_MULA) {
    strcpy(star, ",GPol");
    if ((retflag = swe_fixstar(star, tjd_et, iflag_galequ, x, serr)) == ERR)
      return ERR;
    *daya = swe_degnorm(x[0] - 150 - 6.6666666667);
    return (retflag & SEFLG_EPHMASK);
  }
  if (!(sip->sid_mode & SE_SIDBIT_ECL_DATE)) {

    x[0] = 1;
    x[1] = x[2] = x[3] = x[4] = x[5] = 0;

    if (tjd_et != J2000)
      swi_precess(x, tjd_et, 0, J_TO_J2000);

    t0 = sip->t0;
    if (sip->t0_is_UT)
      t0 += swe_deltat_ex(t0, iflag, serr);
    swi_precess(x, t0, 0, J2000_TO_J);

    eps = swi_epsiln(t0, 0);
    swi_coortrf(x, x, eps);

    swi_cartpol(x, x);

    x[0] = -x[0] * RADTODEG + sip->ayan_t0;
  } else {

    x[0] = swe_degnorm(sip->ayan_t0) * DEGTORAD;
    x[1] = 0; x[2] = 1;

    t0 = sip->t0;
    if (sip->t0_is_UT)
      t0 += swe_deltat_ex(t0, iflag, serr);
    eps = swi_epsiln(t0, 0);

    swi_polcart(x, x);
    swi_coortrf(x, x, -eps);

    if (t0 != J2000)
      swi_precess(x, t0, 0, J_TO_J2000);

    swi_precess(x, tjd_et, 0, J2000_TO_J);

    eps = swi_epsiln(tjd_et, 0);

    swi_coortrf(x, x, eps);
    swi_cartpol(x, x);
    x[0] = swe_degnorm(x[0] * RADTODEG);
  }
  get_aya_correction(iflag, &corr, serr);

  *daya = swe_degnorm(x[0] - corr);

  return iflag;
}

int32 swi_get_ayanamsa_with_speed(double tjd_et, int32 iflag, double *daya, char *serr)
{
  double daya_t2, t2;
  double tintv = 0.001;
  int32 retflag;
  t2 = tjd_et - tintv;
  retflag = swi_get_ayanamsa_ex(t2, iflag, &daya_t2, serr);
  if (retflag == ERR)
    return ERR;
  retflag = swi_get_ayanamsa_ex(tjd_et, iflag, daya, serr);
  if (retflag == ERR)
    return ERR;
  daya[1] = (daya[0] - daya_t2) / tintv;
  return retflag;
}

int32 CALL_CONV swe_get_ayanamsa_ex_ut(double tjd_ut, int32 iflag, double *daya, char *serr)
{
  double deltat;
  int32 retflag = OK;
  int32 epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag == 0) {
    epheflag = SEFLG_SWIEPH;
    iflag |= SEFLG_SWIEPH;
  }
  deltat = swe_deltat_ex(tjd_ut, iflag, serr);

  retflag = swe_get_ayanamsa_ex(tjd_ut + deltat, iflag, daya, serr);

  if ((retflag & SEFLG_EPHMASK) != epheflag) {
    deltat = swe_deltat_ex(tjd_ut, retflag, serr);
    retflag = swe_get_ayanamsa_ex(tjd_ut + deltat, iflag, daya, serr);
  }
  return retflag;
}

double CALL_CONV swe_get_ayanamsa(double tjd_et)
{
  double daya;
  int32 iflag = swi_guess_ephe_flag();

  swi_get_ayanamsa_ex(tjd_et, iflag, &daya, NULL);
  return daya;
}

double CALL_CONV swe_get_ayanamsa_ut(double tjd_ut)
{
  double daya;
  int32 iflag = swi_guess_ephe_flag();
  swi_get_ayanamsa_ex(tjd_ut + swe_deltat_ex(tjd_ut, iflag, NULL), 0, &daya, NULL);
  return daya;
}

int swi_trop_ra2sid_lon(double *xin, double *xout, double *xoutr, int32 iflag)
{
  double x[6], corr;
  int i;
  struct sid_data *sip = &swed.sidd;
  struct epsilon oectmp;
  for (i = 0; i <= 5; i++)
    x[i] = xin[i];
  if (sip->t0 != J2000) {

    swi_precess(x, sip->t0, 0, J2000_TO_J);
    swi_precess(x+3, sip->t0, 0, J2000_TO_J);
  }
  for (i = 0; i <= 5; i++)
    xoutr[i] = x[i];
  calc_epsilon(swed.sidd.t0, iflag, &oectmp);
  swi_coortrf2(x, x, oectmp.seps, oectmp.ceps);
  if (iflag & SEFLG_SPEED)
    swi_coortrf2(x+3, x+3, oectmp.seps, oectmp.ceps);

  swi_cartpol_sp(x, x);

  get_aya_correction(iflag, &corr, NULL);
  x[0] -= sip->ayan_t0 * DEGTORAD;
  x[0] = swe_radnorm(x[0] + corr * DEGTORAD);

  swi_polcart_sp(x, xout);
  return OK;
}

int swi_trop_ra2sid_lon_sosy(double *xin, double *xout, int32 iflag)
{
  double x[6], x0[6], corr;
  int i;
  struct sid_data *sip = &swed.sidd;
  struct epsilon *oe = &swed.oec2000;
  double plane_node = SSY_PLANE_NODE_E2000;
  double plane_incl = SSY_PLANE_INCL;
  for (i = 0; i <= 5; i++)
    x[i] = xin[i];

  swi_coortrf2(x, x, oe->seps, oe->ceps);
  if (iflag & SEFLG_SPEED)
    swi_coortrf2(x+3, x+3, oe->seps, oe->ceps);

  swi_cartpol_sp(x, x);

  x[0] -= plane_node;
  swi_polcart_sp(x, x);
  swi_coortrf(x, x, plane_incl);
  swi_coortrf(x+3, x+3, plane_incl);
  swi_cartpol_sp(x, x);

  x0[0] = 1;
  x0[1] = x0[2] = 0;
  if (sip->t0 != J2000) {

    swi_precess(x0, sip->t0, 0, J_TO_J2000);
  }

  swi_coortrf2(x0, x0, oe->seps, oe->ceps);

  swi_cartpol(x0, x0);

  x0[0] -= plane_node;
  swi_polcart(x0, x0);
  swi_coortrf(x0, x0, plane_incl);
  swi_cartpol(x0, x0);

  x[0] -= x0[0];
  x[0] *= RADTODEG;

  get_aya_correction(iflag, &corr, NULL);
  x[0] -= sip->ayan_t0;
  x[0] = swe_degnorm(x[0] + corr) * DEGTORAD;

  swi_polcart_sp(x, xout);
  return OK;
}

static int app_pos_etc_plan_osc(int ipl, int ipli, int32 iflag, char *serr)
{
  int i, j, niter, retc;
  double xx[6], dx[3], dt, dtsave_for_defl;
  double xearth[6], xsun[6], xmoon[6];
  double xxsv[6], xxsp[3]={0}, xobs[6], xobs2[6];
  double t;
  struct plan_data *pdp = &swed.pldat[ipli];
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed.oec2000;
  int32 epheflag = SEFLG_DEFAULTEPH;
  dt = dtsave_for_defl = 0;
  if (iflag & SEFLG_MOSEPH) {
    epheflag = SEFLG_MOSEPH;
  } else if (iflag & SEFLG_SWIEPH) {
    epheflag = SEFLG_SWIEPH;
  } else if (iflag & SEFLG_JPLEPH) {
    epheflag = SEFLG_JPLEPH;
  }

  for (i = 0; i <= 5; i++)
    xx[i] = pdp->x[i];

  if (iflag & SEFLG_TOPOCTR) {
    if (swed.topd.teval != pedp->teval
      || swed.topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr) != OK)
        return ERR;
    } else {
      for (i = 0; i <= 5; i++)
        xobs[i] = swed.topd.xobs[i];
    }

    for (i = 0; i <= 5; i++)
      xobs[i] = xobs[i] + pedp->x[i];
  } else if (iflag & SEFLG_BARYCTR) {
    for (i = 0; i <= 5; i++)
      xobs[i] = 0;
  } else if (iflag & SEFLG_HELCTR) {
    if (iflag & SEFLG_MOSEPH) {
      for (i = 0; i <= 5; i++)
        xobs[i] = 0;
    } else {
      for (i = 0; i <= 5; i++)
        xobs[i] = psdp->x[i];
    }
  } else {
    for (i = 0; i <= 5; i++)
      xobs[i] = pedp->x[i];
  }

  if (!(iflag & SEFLG_TRUEPOS)) {
    niter = 1;
    if (iflag & SEFLG_SPEED) {

      for (i = 0; i <= 2; i++)
	xxsv[i] = xxsp[i] = xx[i] - xx[i+3];
      for (j = 0; j <= niter; j++) {
	for (i = 0; i <= 2; i++) {
	  dx[i] = xxsp[i];
	  if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR))
	    dx[i] -= (xobs[i] - xobs[i+3]);
	}

	dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
	for (i = 0; i <= 2; i++)
	  xxsp[i] = xxsv[i] - dt * pdp->x[i+3];
      }

      for (i = 0; i <= 2; i++)
	xxsp[i] = xxsv[i] - xxsp[i];
    }

    for (j = 0; j <= niter; j++) {
      for (i = 0; i <= 2; i++) {
	dx[i] = xx[i];
	if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR))
	  dx[i] -= xobs[i];
      }

      dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
      dtsave_for_defl = dt;

      for (i = 0; i <= 2; i++) {
	xx[i] = pdp->x[i] - dt * pdp->x[i+3];
	xx[i+3] = pdp->x[i+3];
      }
    }
    if (iflag & SEFLG_SPEED) {

      for (i = 0; i <= 2; i++)
	xxsp[i] = pdp->x[i] - xx[i] - xxsp[i];
      t = pdp->teval - dt;

      retc = main_planet_bary(t, SEI_EARTH, epheflag, iflag, NO_SAVE, xearth, xearth, xsun, xmoon, serr);
      if (swi_osc_el_plan(t, xx, ipl-SE_FICT_OFFSET, ipli, xearth, xsun, serr) != OK)
	return ERR;
      if (retc != OK)
	return(retc);
      if (iflag & SEFLG_TOPOCTR) {
        if (swi_get_observer(t, iflag | SEFLG_NONUT, NO_SAVE, xobs2, serr) != OK)
          return ERR;
        for (i = 0; i <= 5; i++)
          xobs2[i] += xearth[i];
      } else {
        for (i = 0; i <= 5; i++)
          xobs2[i] = xearth[i];
      }
    }
  }

  for (i = 0; i <= 5; i++)
    xx[i] -= xobs[i];
  if (!(iflag & SEFLG_TRUEPOS)) {

    if (iflag & SEFLG_SPEED)
      for (i = 3; i <= 5; i++)
	xx[i] -= xxsp[i-3];
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOGDEFL))

    swi_deflect_light(xx, dtsave_for_defl, iflag);

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOABERR)) {

    swi_aberr_light(xx, xobs, iflag);

    if (iflag & SEFLG_SPEED)
      for (i = 3; i <= 5; i++)
	xx[i] += xobs[i] - xobs2[i];
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];

  if (!(iflag & SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, pdp->teval, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else
    oe = &swed.oec2000;
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr);
}

void swi_precess_speed(double *xx, double t, int32 iflag, int direction)
{
  struct epsilon *oe;
  double fac, dpre, dpre2;
  double tprec = (t - J2000) / 36525.0;
  int prec_model = swed.astro_models[SE_MODEL_PREC_LONGTERM];
  if (prec_model == 0) prec_model = SEMOD_PREC_DEFAULT;
  if (direction == J2000_TO_J) {
    fac = 1;
    oe = &swed.oec;
  } else {
    fac = -1;
    oe = &swed.oec2000;
  }

  swi_precess(xx+3, t, iflag, direction);

  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
  swi_cartpol_sp(xx, xx);
if (1) {
  if (prec_model == SEMOD_PREC_VONDRAK_2011) {
    swi_ldp_peps(t, &dpre, NULL);
    swi_ldp_peps(t + 1, &dpre2, NULL);
    xx[3] += (dpre2 - dpre) * fac;
  } else {
    xx[3] += (50.290966 + 0.0222226 * tprec) / 3600 / 365.25 * DEGTORAD * fac;

  }
}
  swi_polcart_sp(xx, xx);
  swi_coortrf2(xx, xx, -oe->seps, oe->ceps);
  swi_coortrf2(xx+3, xx+3, -oe->seps, oe->ceps);
}

void swi_nutate(double *xx, int32 iflag, AS_BOOL backward)
{
  int i;
  double x[6], xv[6];
  for (i = 0; i <= 2; i++) {
    if (backward) {
      x[i] = xx[0] * swed.nut.matrix[i][0] +
	     xx[1] * swed.nut.matrix[i][1] +
	     xx[2] * swed.nut.matrix[i][2];
    } else {
      x[i] = xx[0] * swed.nut.matrix[0][i] +
	     xx[1] * swed.nut.matrix[1][i] +
	     xx[2] * swed.nut.matrix[2][i];
    }
  }
  if (iflag & SEFLG_SPEED) {

    for (i = 0; i <= 2; i++) {
      if (backward) {
	x[i+3] = xx[3] * swed.nut.matrix[i][0] +
		 xx[4] * swed.nut.matrix[i][1] +
		 xx[5] * swed.nut.matrix[i][2];
      } else {
	x[i+3] = xx[3] * swed.nut.matrix[0][i] +
		 xx[4] * swed.nut.matrix[1][i] +
		 xx[5] * swed.nut.matrix[2][i];
      }
    }

    for (i = 0; i <= 2; i++) {
      if (backward) {
	xv[i] = xx[0] * swed.nutv.matrix[i][0] +
	       xx[1] * swed.nutv.matrix[i][1] +
	       xx[2] * swed.nutv.matrix[i][2];
      } else {
	xv[i] = xx[0] * swed.nutv.matrix[0][i] +
	       xx[1] * swed.nutv.matrix[1][i] +
	       xx[2] * swed.nutv.matrix[2][i];
      }

      xx[3+i] = x[3+i] + (x[i] - xv[i]) / NUT_SPEED_INTV;
    }
  }

  for (i = 0; i <= 2; i++)
    xx[i] = x[i];
}

static void aberr_light(double *xx, double *xe) {
  int i;
  double xxs[6], v[6], u[6], ru;
  double b_1, f1, f2;
  double v2;
  for (i = 0; i <= 5; i++)
    u[i] = xxs[i] = xx[i];
  ru = sqrt(square_sum(u));
  for (i = 0; i <= 2; i++)
    v[i] = xe[i+3] / 24.0 / 3600.0 / CLIGHT * AUNIT;
  v2 = square_sum(v);
  b_1 = sqrt(1 - v2);
  f1 = dot_prod(u, v) / ru;
  f2 = 1.0 + f1 / (1.0 + b_1);
  for (i = 0; i <= 2; i++)
    xx[i] = (b_1*xx[i] + f2*ru*v[i]) / (1.0 + f1);
}

void swi_aberr_light_ex(double *xx, double *xe, double *xe_dt, double dt, int32 iflag) {
  int i;
  double xxs[6];
  double xx2[6];
  for (i = 0; i <= 5; i++) {
    xxs[i] = xx[i];
  }
  aberr_light(xx, xe);

  if (iflag & SEFLG_SPEED) {
    for (i = 0; i <= 2; i++)
      xx2[i] = xxs[i] - dt * xxs[i + 3];
    aberr_light(xx2, xe_dt);
    for (i = 0; i <= 2; i++) {
      xx[i+3] = (xx[i] - xx2[i]) / dt;
    }
  }
}

void swi_aberr_light(double *xx, double *xe, int32 iflag) {
  int i;
  double xxs[6], v[6], u[6], ru;
  double xx2[6], dx1, dx2;
  double b_1, f1, f2;
  double v2;
  double intv = PLAN_SPEED_INTV;
  for (i = 0; i <= 5; i++)
    u[i] = xxs[i] = xx[i];
  ru = sqrt(square_sum(u));
  for (i = 0; i <= 2; i++)
    v[i] = xe[i+3] / 24.0 / 3600.0 / CLIGHT * AUNIT;
  v2 = square_sum(v);
  b_1 = sqrt(1 - v2);
  f1 = dot_prod(u, v) / ru;
  f2 = 1.0 + f1 / (1.0 + b_1);
  for (i = 0; i <= 2; i++)
    xx[i] = (b_1*xx[i] + f2*ru*v[i]) / (1.0 + f1);
  if (iflag & SEFLG_SPEED) {

    for (i = 0; i <= 2; i++)
      u[i] = xxs[i] - intv * xxs[i+3];
    ru = sqrt(square_sum(u));
    f1 = dot_prod(u, v) / ru;
    f2 = 1.0 + f1 / (1.0 + b_1);
    for (i = 0; i <= 2; i++)
      xx2[i] = (b_1*u[i] + f2*ru*v[i]) / (1.0 + f1);
    for (i = 0; i <= 2; i++) {
      dx1 = xx[i] - xxs[i];
      dx2 = xx2[i] - u[i];
      dx1 -= dx2;
      xx[i+3] += dx1 / intv;
    }
  }
}

void swi_deflect_light(double *xx, double dt, int32 iflag)
{
  int i;
  double xx2[6];
  double u[6], e[6], q[6], ru, re, rq, uq, ue, qe, g1, g2;
#if 1
  double xx3[6], dx1, dx2, dtsp;
#endif
  double xsun[6], xearth[6];
  double sina, sin_sunr, meff_fact;
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  int32 iephe = pedp->iephe;
  for (i = 0; i <= 5; i++)
    xearth[i] = pedp->x[i];
  if (iflag & SEFLG_TOPOCTR)
    for (i = 0; i <= 5; i++)
      xearth[i] += swed.topd.xobs[i];

  for (i = 0; i <= 2; i++)
    u[i] = xx[i];

  if (iephe == SEFLG_JPLEPH || iephe == SEFLG_SWIEPH) {
    for (i = 0; i <= 2; i++)
      e[i] = xearth[i] - psdp->x[i];
  } else {
    for (i = 0; i <= 2; i++)
      e[i] = xearth[i];
  }

  if (iephe == SEFLG_JPLEPH || iephe == SEFLG_SWIEPH) {
    for (i = 0; i <= 2; i++)

      xsun[i] = psdp->x[i] - dt * psdp->x[i+3];
    for (i = 3; i <= 5; i++)
      xsun[i] = psdp->x[i];
  } else {
    for (i = 0; i <= 5; i++)
      xsun[i] = psdp->x[i];
  }
  for (i = 0; i <= 2; i++)
    q[i] = xx[i] + xearth[i] - xsun[i];
  ru = sqrt(square_sum(u));
  rq = sqrt(square_sum(q));
  re = sqrt(square_sum(e));
  for (i = 0; i <= 2; i++) {
    u[i] /= ru;
    q[i] /= rq;
    e[i] /= re;
  }
  uq = dot_prod(u,q);
  ue = dot_prod(u,e);
  qe = dot_prod(q,e);

  sina = sqrt(1 - ue * ue);
  sin_sunr = SUN_RADIUS / re;
  if (sina < sin_sunr) 	{
    meff_fact = meff(sina / sin_sunr);
  } else {
    meff_fact = 1;
  }
  g1 = 2.0 * HELGRAVCONST * meff_fact / CLIGHT / CLIGHT / AUNIT / re;
  g2 = 1.0 + qe;

  for (i = 0; i <= 2; i++)
    xx2[i] = ru * (u[i] + g1/g2 * (uq * e[i] - ue * q[i]));
  if (iflag & SEFLG_SPEED) {

    dtsp = -DEFL_SPEED_INTV;

    for (i = 0; i <= 2; i++)
      u[i] = xx[i] - dtsp * xx[i+3];

    if (iephe == SEFLG_JPLEPH || iephe == SEFLG_SWIEPH) {
      for (i = 0; i <= 2; i++)
	e[i] = xearth[i] - psdp->x[i] -
	       dtsp * (xearth[i+3] - psdp->x[i+3]);
    } else
      for (i = 0; i <= 2; i++)
	e[i] = xearth[i] - dtsp * xearth[i+3];

    for (i = 0; i <= 2; i++)
      q[i] = u[i] + xearth[i] - xsun[i] -
	     dtsp * (xearth[i+3] - xsun[i+3]);
    ru = sqrt(square_sum(u));
    rq = sqrt(square_sum(q));
    re = sqrt(square_sum(e));
    for (i = 0; i <= 2; i++) {
      u[i] /= ru;
      q[i] /= rq;
      e[i] /= re;
    }
    uq = dot_prod(u,q);
    ue = dot_prod(u,e);
    qe = dot_prod(q,e);
    sina = sqrt(1 - ue * ue);
    sin_sunr = SUN_RADIUS / re;
    if (sina < sin_sunr) {
      meff_fact = meff(sina / sin_sunr);
    } else {
      meff_fact = 1;
    }
    g1 = 2.0 * HELGRAVCONST * meff_fact / CLIGHT / CLIGHT / AUNIT / re;
    g2 = 1.0 + qe;
    for (i = 0; i <= 2; i++)
      xx3[i] = ru * (u[i] + g1/g2 * (uq * e[i] - ue * q[i]));
    for (i = 0; i <= 2; i++) {
      dx1 = xx2[i] - xx[i];
      dx2 = xx3[i] - u[i] * ru;
      dx1 -= dx2;
      xx[i+3] += dx1 / dtsp;
    }
  }

  for (i = 0; i <= 2; i++)
    xx[i] = xx2[i];
}

static int app_pos_etc_sun(int32 iflag, char *serr)
{
  int i, j, niter, retc = OK;
  int32 flg1, flg2;
  double xx[6], xxsv[6], dx[3], dt, t = 0;
  double xearth[6], xsun[6], xobs[6];
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed.oec2000;

  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pedp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pedp->xflgs = iflag;
    pedp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }

  if (iflag & SEFLG_TOPOCTR) {
    if (swed.topd.teval != pedp->teval
      || swed.topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr) != OK)
        return ERR;
    } else {
      for (i = 0; i <= 5; i++)
        xobs[i] = swed.topd.xobs[i];
    }

    for (i = 0; i <= 5; i++)
      xobs[i] = xobs[i] + pedp->x[i];
  } else {

    for (i = 0; i <= 5; i++)
      xobs[i] = pedp->x[i];
  }

  if (pedp->iephe == SEFLG_MOSEPH || (iflag & SEFLG_BARYCTR)) {
    for (i = 0; i <= 5; i++)
      xx[i] = xobs[i];
  } else {
    for (i = 0; i <= 5; i++)
      xx[i] = xobs[i] - psdp->x[i];
  }

  if (!(iflag & SEFLG_TRUEPOS)) {

    if (pedp->iephe == SEFLG_JPLEPH || pedp->iephe == SEFLG_SWIEPH
      || (iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {
      for (i = 0; i <= 5; i++) {
        xearth[i] = xobs[i];
	if (pedp->iephe == SEFLG_MOSEPH)
	  xsun[i] = 0;
	else
	  xsun[i] = psdp->x[i];
      }
      niter = 1;
      for (j = 0; j <= niter; j++) {

	for (i = 0; i <= 2; i++) {
	  dx[i] = xearth[i];
	  if (!(iflag & SEFLG_BARYCTR))
	    dx[i] -= xsun[i];
	}

	dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
	t = pedp->teval - dt;

	switch(pedp->iephe) {

	  case SEFLG_JPLEPH:
	    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR))
	      retc = swi_pleph(t, J_EARTH, J_SBARY, xearth, serr);
	    else
	      retc = swi_pleph(t, J_SUN, J_SBARY, xsun, serr);
	    if (retc != OK) {
	      swi_close_jpl_file();
	      swed.jpl_file_is_open = FALSE;
	      return(retc);
	    }
	    break;
	  case SEFLG_SWIEPH:

	    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR)) {
	      retc = sweplan(t, SEI_EARTH, SEI_FILE_PLANET, iflag, NO_SAVE, xearth, NULL, xsun, NULL, serr);
            } else {
	      retc = sweph(t, SEI_SUNBARY, SEI_FILE_PLANET, iflag, NULL, NO_SAVE, xsun, serr);
	    }
	    break;
	  case SEFLG_MOSEPH:
	    if ((iflag & SEFLG_HELCTR) || (iflag & SEFLG_BARYCTR))
	      retc = swi_moshplan(t, SEI_EARTH, NO_SAVE, xearth, xearth, serr);

	    break;
          default:
	    retc = ERR;
	    break;
	}
	if (retc != OK)
	  return(retc);
      }

      for (i = 0; i <= 5; i++) {
        xx[i] = xearth[i];
	if (!(iflag & SEFLG_BARYCTR))
	  xx[i] -= xsun[i];
      }
    }
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR))
    for (i = 0; i <= 5; i++)
      xx[i] = -xx[i];

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOABERR)) {

    swi_aberr_light(xx, xobs, iflag);
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_ICRS) && swi_get_denum(SEI_SUN, iflag) >= 403) {
    swi_bias(xx, t, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];

  if (!(iflag & SEFLG_J2000)) {
    swi_precess(xx, pedp->teval, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, pedp->teval, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else
    oe = &swed.oec2000;
  return app_pos_rest(pedp, iflag, xx, xxsv, oe, serr);
}

static int app_pos_etc_moon(int32 iflag, char *serr)
{
  int i;
  int32 flg1, flg2;
  double xx[6], xxsv[6], xobs[6], xxm[6], xs[6], xe[6], xobs2[6], dt;
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psdp = &swed.pldat[SEI_SUNBARY];
  struct plan_data *pdp = &swed.pldat[SEI_MOON];
  struct epsilon *oe = &swed.oec;
  double t = 0;
  int32 retc;

  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pdp->xflgs = iflag;
    pdp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }

  for (i = 0; i <= 5; i++) {
    xx[i] = pdp->x[i];
    xxm[i] = xx[i];
  }

  for (i = 0; i <= 5; i++)
	xx[i] += pedp->x[i];

  if (iflag & SEFLG_TOPOCTR) {
    if (swed.topd.teval != pdp->teval
      || swed.topd.teval == 0) {
      if (swi_get_observer(pdp->teval, iflag | SEFLG_NONUT, DO_SAVE, xobs, serr) != OK)
        return ERR;
    } else {
      for (i = 0; i <= 5; i++)
        xobs[i] = swed.topd.xobs[i];
    }
    for (i = 0; i <= 5; i++)
      xxm[i] -= xobs[i];
    for (i = 0; i <= 5; i++)
      xobs[i] += pedp->x[i];
  } else if (iflag & SEFLG_BARYCTR) {
    for (i = 0; i <= 5; i++)
      xobs[i] = 0;
    for (i = 0; i <= 5; i++)
      xxm[i] += pedp->x[i];
  } else if (iflag & SEFLG_HELCTR) {
    for (i = 0; i <= 5; i++)
      xobs[i] = psdp->x[i];
    for (i = 0; i <= 5; i++)
      xxm[i] += pedp->x[i] - psdp->x[i];
  } else {
    for (i = 0; i <= 5; i++)
      xobs[i] = pedp->x[i];
  }

  t = pdp->teval;
  if ((iflag & SEFLG_TRUEPOS) == 0) {
    dt = sqrt(square_sum(xxm)) * AUNIT / CLIGHT / 86400.0;
    t = pdp->teval - dt;
    switch(pdp->iephe) {
      case SEFLG_JPLEPH:
        retc = swi_pleph(t, J_MOON, J_EARTH, xx, serr);
        if (retc == OK)
          retc = swi_pleph(t, J_EARTH, J_SBARY, xe, serr);
        if (retc == OK && (iflag & SEFLG_HELCTR))
          retc = swi_pleph(t, J_SUN, J_SBARY, xs, serr);
        if (retc != OK) {
	      swi_close_jpl_file();
	      swed.jpl_file_is_open = FALSE;
        }
	for (i = 0; i <= 5; i++)
	  xx[i] += xe[i];
	    break;
      case SEFLG_SWIEPH:
        retc = sweplan(t, SEI_MOON, SEI_FILE_MOON, iflag, NO_SAVE, xx, xe, xs, NULL, serr);
        if (retc != OK)
          return(retc);
	for (i = 0; i <= 5; i++)
	  xx[i] += xe[i];
	    break;
      case SEFLG_MOSEPH:

        for (i = 0; i <= 2; i++) {
          xx[i] -= dt * xx[i+3];
          xe[i] = pedp->x[i] - dt * pedp->x[i+3];
		  xe[i+3] = pedp->x[i+3];
	  xs[i] = 0;
	  xs[i+3] = 0;
        }
        break;
    }
    if (iflag & SEFLG_TOPOCTR) {
      if (swi_get_observer(t, iflag | SEFLG_NONUT, NO_SAVE, xobs2, NULL) != OK)
	  return ERR;
      for (i = 0; i <= 5; i++)
	xobs2[i] += xe[i];
    } else if (iflag & SEFLG_BARYCTR) {
      for (i = 0; i <= 5; i++)
	xobs2[i] = 0;
    } else if (iflag & SEFLG_HELCTR) {
      for (i = 0; i <= 5; i++)
	xobs2[i] = xs[i];
    } else {
      for (i = 0; i <= 5; i++)
	xobs2[i] = xe[i];
    }
  }

  for (i = 0; i <= 5; i++)
    xx[i] -= xobs[i];

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOABERR)) {

    swi_aberr_light(xx, xobs, iflag);

#if 1
    if (iflag & SEFLG_SPEED)
      for (i = 3; i <= 5; i++)
        xx[i] += xobs[i] - xobs2[i];
#endif
  }

  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_ICRS) && swi_get_denum(SEI_MOON, iflag) >= 403) {
    swi_bias(xx, t, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];

  if (!(iflag & SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, pdp->teval, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else
    oe = &swed.oec2000;
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr);
}

static int app_pos_etc_sbar(int32 iflag, char *serr)
{
  int i;
  double xx[6], xxsv[6], dt;
  struct plan_data *psdp = &swed.pldat[SEI_EARTH];
  struct plan_data *psbdp = &swed.pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed.oec;

  for (i = 0; i <= 5; i++)
    xx[i] = psbdp->x[i];

  if (!(iflag & SEFLG_TRUEPOS)) {
    dt = sqrt(square_sum(xx)) * AUNIT / CLIGHT / 86400.0;
    for (i = 0; i <= 2; i++)
      xx[i] -= dt * xx[i+3];
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_ICRS) && swi_get_denum(SEI_SUN, iflag) >= 403) {
    swi_bias(xx, psdp->teval, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];

  if (!(iflag & SEFLG_J2000)) {
    swi_precess(xx, psbdp->teval, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, psbdp->teval, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else
    oe = &swed.oec2000;
  return app_pos_rest(psdp, iflag, xx, xxsv, oe, serr);
}

static int app_pos_etc_mean(int ipl, int32 iflag, char *serr)
{
  int i;
  int32 flg1, flg2;
  double xx[6], xxsv[6];
  struct plan_data *pdp = &swed.nddat[ipl];
  struct epsilon *oe;

  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pdp->xflgs = iflag;
    pdp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }
  for (i = 0; i <= 5; i++)
    xx[i] = pdp->x[i];

  swi_polcart_sp(xx, xx);
  swi_coortrf2(xx, xx, -swed.oec.seps, swed.oec.ceps);
  swi_coortrf2(xx+3, xx+3, -swed.oec.seps, swed.oec.ceps);
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (((iflag & SEFLG_SIDEREAL)
    && (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0))
      || (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE)) {
    for (i = 0; i <= 5; i++)
      xxsv[i] = xx[i];

    if (pdp->teval != J2000) {
      swi_precess(xxsv, pdp->teval, iflag, J_TO_J2000);
      if (iflag & SEFLG_SPEED)
        swi_precess_speed(xxsv, pdp->teval, iflag, J_TO_J2000);
    }
  }

  if (iflag & SEFLG_J2000) {
    swi_precess(xx, pdp->teval, iflag, J_TO_J2000);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, pdp->teval, iflag, J_TO_J2000);
    oe = &swed.oec2000;
  } else
    oe = &swed.oec;
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr);
}

static int get_new_segment(double tjd, int ipli, int ifno, char *serr)
{
  int i, j, k, m, n, o, icoord, retc;
  int32 iseg;
  int32 fpos;
  int nsizes, nsize[6];
  int nco;
  int idbl;
  unsigned char c[4];
  struct plan_data *pdp = &swed.pldat[ipli];
  struct file_data *fdp = &swed.fidat[ifno];
  FILE *fp = fdp->fptr;
  int freord  = (int) fdp->iflg & SEI_FILE_REORD;
  int fendian = (int) fdp->iflg & SEI_FILE_LITENDIAN;
  uint32 longs[MAXORD+1];

  iseg = (int32) ((tjd - pdp->tfstart) / pdp->dseg);

  pdp->tseg0 = pdp->tfstart + iseg * pdp->dseg;
  pdp->tseg1 = pdp->tseg0 + pdp->dseg;

  fpos = pdp->lndx0 + iseg * 3;
  retc = do_fread((void *) &fpos, 3, 1, 4, fp, fpos, freord, fendian, ifno, serr);
  if (retc != OK)
    goto return_error_gns;
  fseek(fp, fpos, SEEK_SET);

  if (pdp->segp == NULL)
    pdp->segp = (double *) malloc((size_t) pdp->ncoe * 3 * 8);
  memset((void *) pdp->segp, 0, (size_t) pdp->ncoe * 3 * 8);

  for (icoord = 0; icoord < 3; icoord++) {
    idbl = icoord * pdp->ncoe;

    retc = do_fread((void *) &c[0], 1, 2, 1, fp, SEI_CURR_FPOS, freord, fendian, ifno, serr);
    if (retc != OK)
      goto return_error_gns;
    if (c[0] & 128) {
      nsizes = 6;
      retc = do_fread((void *) (c+2), 1, 2, 1, fp, SEI_CURR_FPOS, freord, fendian, ifno, serr);
      if (retc != OK)
	goto return_error_gns;
      nsize[0] = (int) c[1] / 16;
      nsize[1] = (int) c[1] % 16;
      nsize[2] = (int) c[2] / 16;
      nsize[3] = (int) c[2] % 16;
      nsize[4] = (int) c[3] / 16;
      nsize[5] = (int) c[3] % 16;
      nco = nsize[0] + nsize[1] + nsize[2] + nsize[3] + nsize[4] + nsize[5];
    } else {
      nsizes = 4;
      nsize[0] = (int) c[0] / 16;
      nsize[1] = (int) c[0] % 16;
      nsize[2] = (int) c[1] / 16;
      nsize[3] = (int) c[1] % 16;
      nco = nsize[0] + nsize[1] + nsize[2] + nsize[3];
    }

    if (nco > pdp->ncoe) {
      if (serr != NULL) {
	sprintf(serr, "error in ephemeris file: %d coefficients instead of %d. ", nco, pdp->ncoe);
	if (strlen(serr) + strlen(fdp->fnam) < AS_MAXCH - 1) {
	  sprintf(serr, "error in ephemeris file %s: %d coefficients instead of %d. ", fdp->fnam, nco, pdp->ncoe);
	}
      }
      free(pdp->segp);
      pdp->segp = NULL;
      return (ERR);
    }

    for (i = 0; i < nsizes; i++) {
      if (nsize[i] == 0)
	continue;
      if (i < 4) {
	j = (4 - i);
	k = nsize[i];
	retc = do_fread((void *) &longs[0], j, k, 4, fp, SEI_CURR_FPOS, freord, fendian, ifno, serr);
	if (retc != OK)
	  goto return_error_gns;
	for (m = 0; m < k; m++, idbl++) {
	  if (longs[m] & 1)
	    pdp->segp[idbl] = -(((longs[m]+1) / 2) / 1e+9 * pdp->rmax / 2);
	  else
	    pdp->segp[idbl] = (longs[m] / 2) / 1e+9 * pdp->rmax / 2;
	}
      } else if (i == 4) {
	j = 1;
	k = (nsize[i] + 1) / 2;
	retc = do_fread((void *) longs, j, k, 4, fp, SEI_CURR_FPOS, freord, fendian, ifno, serr);
	if (retc != OK)
	  goto return_error_gns;
	for (m = 0, j = 0;
	     m < k && j < nsize[i];
	     m++) {
	  for (n = 0, o = 16;
	       n < 2 && j < nsize[i];
	       n++, j++, idbl++, longs[m] %= o, o /= 16) {
	    if (longs[m] & o)
	      pdp->segp[idbl] =
		   -(((longs[m]+o) / o / 2) * pdp->rmax / 2 / 1e+9);
	    else
	      pdp->segp[idbl] = (longs[m] / o / 2) * pdp->rmax / 2 / 1e+9;
	  }
	}
      } else if (i == 5) {
	j = 1;
	k = (nsize[i] + 3) / 4;
	retc = do_fread((void *) longs, j, k, 4, fp, SEI_CURR_FPOS, freord, fendian, ifno, serr);
	if (retc != OK)
	  goto return_error_gns;
	for (m = 0, j = 0;
	     m < k && j < nsize[i];
	     m++) {
	  for (n = 0, o = 64;
	       n < 4 && j < nsize[i];
	       n++, j++, idbl++, longs[m] %= o, o /= 4) {
	    if (longs[m] & o)
	      pdp->segp[idbl] =
		   -(((longs[m]+o) / o / 2) * pdp->rmax / 2 / 1e+9);
	    else
	      pdp->segp[idbl] = (longs[m] / o / 2) * pdp->rmax / 2 / 1e+9;
	  }
	}
      }
    }
  }
  return(OK);
return_error_gns:
  fclose(fdp->fptr);

  fdp->fptr = NULL;
  free_planets();
  return ERR;
}

static int read_const(int ifno, char *serr)
{
  char *c, c2, *sp;
  char s[AS_MAXCH*2], s2[AS_MAXCH];
  char sastnam[41];
  int i, ipli, kpl;
  int retc;
  int fendian, freord;
  int lastnam = 19;
  FILE *fp;
  int32 lng;
  uint32 ulng;
  int32 flen, fpos;
  short nplan;
  int32 testendian;
  double doubles[20];
  struct plan_data *pdp;
  struct file_data *fdp = &swed.fidat[ifno];
  char *serr_file_damage = "Ephemeris file %s is damaged (0%s). ";
  char *smsg = "";
  int nbytes_ipl = 2;
  fp = fdp->fptr;

  sp = fgets(s, AS_MAXCH, fp);
  if (sp == NULL || strstr(sp, "\r\n") == NULL) {
    goto file_damage;
  }
  sp = strchr(s, '\r');
  *sp = '\0';
  sp = s;
  while (isdigit((int) *sp) == 0 && *sp != '\0')
    sp++;
  if (*sp == '\0') {
    smsg = "a";
    goto file_damage;
  }

  fdp->fversion = atoi(sp);

  sp = fgets(s, AS_MAXCH, fp);
  if (sp == NULL || strstr(sp, "\r\n") == NULL) {
    smsg = "b";
    goto file_damage;
  }

  sp = strrchr(fdp->fnam, (int) *DIR_GLUE);
  if (sp == NULL) {
    sp = fdp->fnam;
  } else {
    sp++;
  }
  strcpy(s2, sp);

  for (sp = s2; *sp != '\0'; sp++)
    *sp = tolower((int) *sp);

  sp = s + strlen(s) - 1;
  while (*sp == '\n' || *sp == '\r' || *sp == ' ') {
    *sp = '\0';
    sp--;
  }
  for (sp = s; *sp != '\0'; sp++)
    *sp = tolower((int) *sp);
  if (strcmp(s2, s) != 0) {
    if (serr != NULL) {
      sprintf(serr, "Ephemeris file name '%s' wrong; rename '%s' ", s2, s);
    }
    goto return_error;
  }

  sp = fgets(s, AS_MAXCH, fp);
  if (sp == NULL || strstr(sp, "\r\n") == NULL) {
    smsg = "c";
    goto file_damage;
  }

  if (ifno == SEI_FILE_ANY_AST) {
    sp = fgets(s, AS_MAXCH * 2, fp);
    if (sp == NULL || strstr(sp, "\r\n") == NULL) {
      smsg = "d";
      goto file_damage;
    }

    while(*sp == ' ') sp++;
    while(isdigit((int) *sp)) sp++;
    sp++;
    i = (int) (sp - s);
    strncpy(sastnam, s, lastnam+i);
    *(sastnam+lastnam+i) = '\0';

    strcpy(swed.astelem, s);

    swed.ast_H = atof(s + 35 + i);
    swed.ast_G = atof(s + 42 + i);
    if (swed.ast_G == 0) swed.ast_G = 0.15;

    strncpy(s2, s+51+i, 7);
    *(s2 + 7) = '\0';
    swed.ast_diam = atof(s2);
    if (swed.ast_diam == 0) {

      swed.ast_diam = 1329/sqrt(0.15) * pow(10, -0.2 * swed.ast_H);
    }
  }

  if (fread((void *) &testendian, 4, 1, fp) != 1) {
    smsg = "e";
    goto file_damage;
  }

  if (testendian == SEI_FILE_TEST_ENDIAN) {
    freord = SEI_FILE_NOREORD;
  } else {
    freord = SEI_FILE_REORD;
    sp = (char *) &lng;
    c = (char *) &testendian;
    for (i = 0; i < 4; i++)
      *(sp+i) = *(c+3-i);
    if (lng != SEI_FILE_TEST_ENDIAN) {
      smsg = "f";
      goto file_damage;
    }
  }

  c = (char *) &testendian;
  c2 = SEI_FILE_TEST_ENDIAN / 16777216L;
  if (*c == c2) {
    fendian = SEI_FILE_BIGENDIAN;
  } else {
    fendian = SEI_FILE_LITENDIAN;
  }
  fdp->iflg = (int32) freord | fendian;

  retc = do_fread((void *) &lng, 4, 1, 4, fp, SEI_CURR_FPOS, freord,
fendian, ifno, serr);
  if (retc != OK)
    goto return_error;
  fpos = (int32) ftell(fp);
  if (fseek(fp, 0L, SEEK_END) != 0) {
    smsg = "g";
    goto file_damage;
  }
  flen = (int32) ftell(fp);
  if (lng != flen) {
    smsg = "h";
    goto file_damage;
  }

  retc = do_fread((void *) &fdp->sweph_denum, 4, 1, 4, fp, fpos, freord,
fendian, ifno, serr);
  if (retc != OK)
    goto return_error;

  retc = do_fread((void *) &fdp->tfstart, 8, 1, 8, fp, SEI_CURR_FPOS,
freord, fendian, ifno, serr);
  if (retc != OK)
    goto return_error;
  retc = do_fread((void *) &fdp->tfend, 8, 1, 8, fp, SEI_CURR_FPOS, freord,
fendian, ifno, serr);
  if (retc != OK)
    goto return_error;

  retc = do_fread((void *) &nplan, 2, 1, 2, fp, SEI_CURR_FPOS, freord, fendian, ifno, serr);
  if (retc != OK)
    goto return_error;
  if (nplan > 256) {
    nbytes_ipl = 4;
    nplan %= 256;
  }
  if (nplan < 1 || nplan > 20) {
    smsg = "i";
    goto file_damage;
  }
  fdp->npl = nplan;

  retc = do_fread((void *) fdp->ipl, nbytes_ipl, (int) nplan, sizeof(int), fp, SEI_CURR_FPOS,
freord, fendian, ifno, serr);
  if (retc != OK)
    goto return_error;

  if (ifno == SEI_FILE_ANY_AST) {
    char sastno[12];
    int j;

    j = 4;
    while (sastnam[j] != ' ' && j < 10)
      j++;
    strncpy(sastno, sastnam, j);
    sastno[j] = '\0';
    i = (int) atol(sastno);
    if (i == fdp->ipl[0] - SE_AST_OFFSET ||
        i == fdp->ipl[0]
	) {

      strncpy(fdp->astnam, sastnam+j+1, lastnam);
      fdp->astnam[lastnam] = '\0';

      if (fread((void *) s, 30, 1, fp) != 1) {
	smsg = "j";
        goto file_damage;
      }
    } else {

      if (fread((void *) fdp->astnam, 30, 1, fp) != 1) {
	smsg = "k";
        goto file_damage;
      }
    }

    i = (int) (strlen(fdp->astnam) - 1);
    if (i < 0)
      i = 0;
    sp = fdp->astnam + i;
    while(*sp == ' ') {
      sp--;
    }
    sp[1] = '\0';
    if ((sp = strstr(fdp->astnam, "  ")) != NULL)
      *sp = '\0';
  }

  fpos = (int32) ftell(fp);

  retc = do_fread((void *) &ulng, 4, 1, 4, fp, SEI_CURR_FPOS, freord,
fendian, ifno, serr);
  if (retc != OK)
    goto return_error;

  fseek(fp, 0L, SEEK_SET);

  if (fpos - 1 > 2 * AS_MAXCH) {
    smsg = "l";
    goto file_damage;
  }
  if (fread((void *) s, (size_t) fpos, 1, fp) != 1) {
    smsg = "m";
    goto file_damage;
  }
#if 1
  if (swi_crc32((unsigned char *) s, (int) fpos) != ulng) {
    smsg = "n";
    goto file_damage;
  }
#endif
  fseek(fp, fpos+4, SEEK_SET);

  retc = do_fread((void *) &doubles[0], 8, 5, 8, fp, SEI_CURR_FPOS, freord,
fendian, ifno, serr);
  if (retc != OK)
    goto return_error;
  swed.gcdat.clight       = doubles[0];
  swed.gcdat.aunit        = doubles[1];
  swed.gcdat.helgravconst = doubles[2];
  swed.gcdat.ratme        = doubles[3];
  swed.gcdat.sunradius    = doubles[4];

  for (kpl = 0; kpl < fdp->npl; kpl++) {

    ipli = fdp->ipl[kpl];
    if (ipli >= SE_AST_OFFSET) {
      pdp = &swed.pldat[SEI_ANYBODY];
    } else if (ipli >= SE_PLMOON_OFFSET) {
      pdp = &swed.pldat[SEI_ANYBODY];
    } else {
      pdp = &swed.pldat[ipli];
    }
    pdp->ibdy = ipli;

    retc = do_fread((void *) &pdp->lndx0, 4, 1, 4, fp, SEI_CURR_FPOS,
freord, fendian, ifno, serr);
    if (retc != OK)
      goto return_error;

    retc = do_fread((void *) &pdp->iflg, 1, 1, sizeof(int32), fp,
SEI_CURR_FPOS, freord, fendian, ifno, serr);
    if (retc != OK)
      goto return_error;

    retc = do_fread((void *) &pdp->ncoe, 1, 1, sizeof(int), fp,
SEI_CURR_FPOS, freord, fendian, ifno, serr);
    if (retc != OK)
      goto return_error;

    retc = do_fread((void *) &lng, 4, 1, 4, fp, SEI_CURR_FPOS, freord,
fendian, ifno, serr);
    if (retc != OK)
      goto return_error;
    pdp->rmax = lng / 1000.0;

    if (ipli >= SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET) {
      if ((ipli % 100) == 99 || (ipli - 9000) / 100 == SE_MARS)
	pdp->rmax = lng / 1000000.0;
    }

    retc = do_fread((void *) doubles, 8, 10, 8, fp, SEI_CURR_FPOS, freord,
fendian, ifno, serr);
    if (retc != OK)
      goto return_error;
    pdp->tfstart  = doubles[0];
    pdp->tfend    = doubles[1];
    pdp->dseg     = doubles[2];
    pdp->nndx     = (int32) ((doubles[1] - doubles[0] + 0.1) /doubles[2]);
    pdp->telem    = doubles[3];
    pdp->prot     = doubles[4];
    pdp->dprot    = doubles[5];
    pdp->qrot     = doubles[6];
    pdp->dqrot    = doubles[7];
    pdp->peri     = doubles[8];
    pdp->dperi    = doubles[9];

    if (pdp->iflg & SEI_FLG_ELLIPSE) {
      if (pdp->refep != NULL) {
        free((void *) pdp->refep);
	pdp->refep = NULL;
        if (pdp->segp != NULL) {
          free((void *) pdp->segp);
          pdp->segp = NULL;
        }
      }
      pdp->refep = (double *) malloc((size_t) pdp->ncoe * 2 * 8);
      retc = do_fread((void *) pdp->refep, 8, 2*pdp->ncoe, 8, fp,
SEI_CURR_FPOS, freord, fendian, ifno, serr);
      if (retc != OK) {
	free(pdp->refep);
	pdp->refep = NULL;
	goto return_error;
      }
    }
  }
  return(OK);
file_damage:
  if (serr != NULL) {
    *serr = '\0';
    if (strlen(serr_file_damage) + strlen(fdp->fnam) + strlen(smsg) < AS_MAXCH) {
      sprintf(serr, serr_file_damage, fdp->fnam, smsg);
    }
  }
return_error:
  fclose(fdp->fptr);

  fdp->fptr = NULL;
  free_planets();
  return(ERR);
}

static int do_fread(void *trg, int size, int count, int corrsize, FILE *fp, int32 fpos, int freord, int fendian, int ifno, char *serr)
{
  int i, j, k;
  int totsize;
  unsigned char space[1000];
  unsigned char *targ = (unsigned char *) trg;
  totsize = size * count;
  if (fpos >= 0)
    fseek(fp, fpos, SEEK_SET);

  if (!freord && size == corrsize) {
    if (fread((void *) targ, (size_t) totsize, 1, fp) == 0) {
      if (serr != NULL) {
	strcpy(serr, "Ephemeris file is damaged (1). ");
	if (strlen(serr) + strlen(swed.fidat[ifno].fnam) < AS_MAXCH - 1) {
	  sprintf(serr, "Ephemeris file %s is damaged (2).", swed.fidat[ifno].fnam);
	}
      }
      return(ERR);
    } else
      return(OK);
  } else {
    if (fread((void *) &space[0], (size_t) totsize, 1, fp) == 0) {
      if (serr != NULL) {
	strcpy(serr, "Ephemeris file is damaged (3). ");
	if (strlen(serr) + strlen(swed.fidat[ifno].fnam) < AS_MAXCH - 1) {
	  sprintf(serr, "Ephemeris file %s is damaged (4).", swed.fidat[ifno].fnam);
	}
      }
      return(ERR);
    }
    if (size != corrsize) {
      memset((void *) targ, 0, (size_t) count * corrsize);
    }
    for(i = 0; i < count; i++) {
      for (j = size-1; j >= 0; j--) {
	if (freord) {
	  k = size-j-1;
	} else {
	  k = j;
	}
        if (size != corrsize) {
          if ((fendian == SEI_FILE_BIGENDIAN && !freord) ||
              (fendian == SEI_FILE_LITENDIAN &&  freord))
	    k += corrsize - size;
	}
        targ[i*corrsize+k] = space[i*size+j];
      }
    }
  }
  return(OK);
}

static void rot_back(int ipli)
{
  int i;
  double t, tdiff;
  double qav, pav, dn;
  double omtild, com, som, cosih2;
  double x[MAXORD+1][3];
  double uix[3], uiy[3], uiz[3];
  double xrot, yrot, zrot;
  double *chcfx, *chcfy, *chcfz;
  double *refepx, *refepy;

  double seps2000 = 0.39777715572793088;
  double ceps2000 = 0.91748206215761929;
  struct plan_data *pdp = &swed.pldat[ipli];
  int nco = pdp->ncoe;
  t = pdp->tseg0 + pdp->dseg / 2;
  chcfx = pdp->segp;
  chcfy = chcfx + nco;
  chcfz = chcfx + 2 * nco;
  tdiff= (t - pdp->telem) / 365250.0;
  if (ipli == SEI_MOON) {
    dn = pdp->prot + tdiff * pdp->dprot;
    i = (int) (dn / TWOPI);
    dn -= i * TWOPI;
    qav = (pdp->qrot + tdiff * pdp->dqrot) * cos(dn);
    pav = (pdp->qrot + tdiff * pdp->dqrot) * sin(dn);
  } else {
    qav = pdp->qrot + tdiff * pdp->dqrot;
    pav = pdp->prot + tdiff * pdp->dprot;
  }

  for (i = 0; i < nco; i++) {
    x[i][0] = chcfx[i];
    x[i][1] = chcfy[i];
    x[i][2] = chcfz[i];
  }
  if (pdp->iflg & SEI_FLG_ELLIPSE) {
    refepx = pdp->refep;
    refepy = refepx + nco;
    omtild = pdp->peri + tdiff * pdp->dperi;
    i = (int) (omtild / TWOPI);
    omtild -= i * TWOPI;
    com = cos(omtild);
    som = sin(omtild);

    for (i = 0; i < nco; i++) {
      x[i][0] = chcfx[i] + com * refepx[i] - som * refepy[i];
      x[i][1] = chcfy[i] + com * refepy[i] + som * refepx[i];
    }
  }

  cosih2 = 1.0 / (1.0 + qav * qav + pav * pav);

  uiz[0] = 2.0 * pav * cosih2;
  uiz[1] = -2.0 * qav * cosih2;
  uiz[2] = (1.0 - qav * qav - pav * pav) * cosih2;

  uix[0] = (1.0 + qav * qav - pav * pav) * cosih2;
  uix[1] = 2.0 * qav * pav * cosih2;
  uix[2] = -2.0 * pav * cosih2;

  uiy[0] =2.0 * qav * pav * cosih2;
  uiy[1] =(1.0 - qav * qav + pav * pav) * cosih2;
  uiy[2] =2.0 * qav * cosih2;

  for (i = 0; i < nco; i++) {
    xrot = x[i][0] * uix[0] + x[i][1] * uiy[0] + x[i][2] * uiz[0];
    yrot = x[i][0] * uix[1] + x[i][1] * uiy[1] + x[i][2] * uiz[1];
    zrot = x[i][0] * uix[2] + x[i][1] * uiy[2] + x[i][2] * uiz[2];
    if (fabs(xrot) + fabs(yrot) + fabs(zrot) >= 1e-14)
      pdp->neval = i;
    x[i][0] = xrot;
    x[i][1] = yrot;
    x[i][2] = zrot;
    if (ipli == SEI_MOON) {

      x[i][1] = ceps2000 * yrot - seps2000 * zrot;
      x[i][2] = seps2000 * yrot + ceps2000 * zrot;
    }
  }
  for (i = 0; i < nco; i++) {
    chcfx[i] = x[i][0];
    chcfy[i] = x[i][1];
    chcfz[i] = x[i][2];
  }
}

static void embofs(double *xemb, double *xmoon)
{
  int i;
  for (i = 0; i <= 2; i++)
    xemb[i] -= xmoon[i] / (EARTH_MOON_MRAT + 1.0);
}

static void nut_matrix(struct nut *nu, struct epsilon *oe)
{
  double psi, eps;
  double sinpsi, cospsi, sineps, coseps, sineps0, coseps0;
  psi = nu->nutlo[0];
  eps = oe->eps + nu->nutlo[1];
  sinpsi = sin(psi);
  cospsi = cos(psi);
  sineps0 = oe->seps;
  coseps0 = oe->ceps;
  sineps = sin(eps);
  coseps = cos(eps);
  nu->matrix[0][0] = cospsi;
  nu->matrix[0][1] = sinpsi * coseps;
  nu->matrix[0][2] = sinpsi * sineps;
  nu->matrix[1][0] = -sinpsi * coseps0;
  nu->matrix[1][1] = cospsi * coseps * coseps0 + sineps * sineps0;
  nu->matrix[1][2] = cospsi * sineps * coseps0 - coseps * sineps0;
  nu->matrix[2][0] = -sinpsi * sineps0;
  nu->matrix[2][1] = cospsi * coseps * sineps0 - sineps * coseps0;
  nu->matrix[2][2] = cospsi * sineps * sineps0 + coseps * coseps0;
}

static int lunar_osc_elem(double tjd, int ipl, int32 iflag, char *serr)
{
  int i, j, istart;
  int ipli = SEI_MOON;
  int32 epheflag = SEFLG_DEFAULTEPH;
  int retc = ERR;
  int32 flg1, flg2;
  double daya[2];
  struct plan_data *ndp, *ndnp, *ndap;
  struct epsilon *oe;
  double speed_intv = NODE_CALC_INTV;
  double a, b;
  double xpos[3][6], xx[3][6], xxa[3][6], xnorm[6], r[6];
  double *xp;
  double rxy, rxyz, t, dt, fac, sgn;
  double sinnode, cosnode, sinincl, cosincl, sinu, cosu, sinE, cosE;
  double uu, ny, sema, ecce, Gmsm, c2, v2, pp;
  int32 speedf1, speedf2;
#ifdef SID_TNODE_FROM_ECL_T0
  struct sid_data *sip = &swed.sidd;
  struct epsilon oectmp;
#endif
  oe = &swed.oec;
#ifdef SID_TNODE_FROM_ECL_T0
  if (iflag & SEFLG_SIDEREAL) {
    calc_epsilon(sip->t0, iflag, &oectmp);
    oe = &oectmp;
  } else if (iflag & SEFLG_J2000) {
    oe = &swed.oec2000;
  }
#endif
  ndp = &swed.nddat[ipl];

  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = ndp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  speedf1 = ndp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == ndp->teval
	&& tjd != 0
	&& flg1 == flg2
	&& (!speedf2 || speedf1)) {
    ndp->xflgs = iflag;
    ndp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }

  if (iflag & SEFLG_MOSEPH) {
    epheflag = SEFLG_MOSEPH;
  } else if (iflag & SEFLG_SWIEPH) {
    epheflag = SEFLG_SWIEPH;
  } else if (iflag & SEFLG_JPLEPH) {
    epheflag = SEFLG_JPLEPH;
  }

  swed.pldat[SEI_MOON].teval = 0;
  if (iflag & SEFLG_SPEED) {
    istart = 0;
  } else {
    istart = 2;
  }
  if (serr != NULL)
    *serr = '\0';
  three_positions:
  switch(epheflag) {
    case SEFLG_JPLEPH:
      speed_intv = NODE_CALC_INTV;
      for (i = istart; i <= 2; i++) {
	if (i == 0) {
	  t = tjd - speed_intv;
        } else if (i == 1) {
	  t = tjd + speed_intv;
        } else  {
	  t = tjd;
	}
	xp = xpos[i];
	retc = jplplan(t, ipli, iflag, NO_SAVE, xp, NULL, NULL, serr);

	if (retc == ERR)
	  return(ERR);

	if ((iflag & SEFLG_TRUEPOS) == 0 && retc >= OK) {
	  dt = sqrt(square_sum(xpos[i])) * AUNIT / CLIGHT / 86400.0;
	  retc = jplplan(t-dt, ipli, iflag, NO_SAVE, xpos[i], NULL, NULL, serr);

	  if (retc == ERR)
	    return(ERR);
        }

	if (retc == NOT_AVAILABLE) {
	  iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_SWIEPH;
	  epheflag = SEFLG_SWIEPH;
	  if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	    strcat(serr, " \ntrying Swiss Eph; ");
	  break;
	} else if (retc == BEYOND_EPH_LIMITS) {
	  if (tjd > MOSHLUEPH_START && tjd < MOSHLUEPH_END) {
	    iflag = (iflag & ~SEFLG_JPLEPH) | SEFLG_MOSEPH;
	    epheflag = SEFLG_MOSEPH;
	    if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	      strcat(serr, " \nusing Moshier Eph; ");
	    break;
	  } else
	    return ERR;
	}

	retc = swi_plan_for_osc_elem(iflag|SEFLG_SPEED, t, xpos[i]);
      }
      break;
    case SEFLG_SWIEPH:
      speed_intv = NODE_CALC_INTV;
      for (i = istart; i <= 2; i++) {
	if (i == 0) {
	  t = tjd - speed_intv;
        } else if (i == 1) {
	  t = tjd + speed_intv;
        } else  {
	  t = tjd;
	}
	retc = swemoon(t, iflag | SEFLG_SPEED, NO_SAVE, xpos[i], serr);
	if (retc == ERR)
	  return(ERR);

	if ((iflag & SEFLG_TRUEPOS) == 0 && retc >= OK) {
	  dt = sqrt(square_sum(xpos[i])) * AUNIT / CLIGHT / 86400.0;
	  retc = swemoon(t-dt, iflag | SEFLG_SPEED, NO_SAVE, xpos[i], serr);
	  if (retc == ERR)
	    return(ERR);
        }
	if (retc == NOT_AVAILABLE) {
	  if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
	    iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
	    epheflag = SEFLG_MOSEPH;
	    if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
	      strcat(serr, " \nusing Moshier eph.; ");
	    break;
	  } else
	    return ERR;
	}

	retc = swi_plan_for_osc_elem(iflag|SEFLG_SPEED, t, xpos[i]);
      }
      break;
    case SEFLG_MOSEPH:

      speed_intv = NODE_CALC_INTV_MOSH;
      for (i = istart; i <= 2; i++) {
	if (i == 0) {
	  t = tjd - speed_intv;
        } else if (i == 1) {
	  t = tjd + speed_intv;
        } else  {
	  t = tjd;
	}
	retc = swi_moshmoon(t, NO_SAVE, xpos[i], serr);
	if (retc == ERR)
	  return(retc);

	retc = swi_plan_for_osc_elem(iflag|SEFLG_SPEED, t, xpos[i]);
      }
      break;
    default:
      break;
  }
  if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS)
    goto three_positions;

  ndnp = &swed.nddat[SEI_TRUE_NODE];

  for (i = istart; i <= 2; i++) {
    if (fabs(xpos[i][5]) < 1e-15)
      xpos[i][5] = 1e-15;
    fac = xpos[i][2] / xpos[i][5];
    sgn = xpos[i][5] / fabs(xpos[i][5]);
    for (j = 0; j <= 2; j++)
      xx[i][j] = (xpos[i][j] - fac * xpos[i][j+3]) * sgn;
  }

  for (i = 0; i <= 2; i++) {
    ndnp->x[i] = xx[2][i];
    if (iflag & SEFLG_SPEED) {
      b = (xx[1][i] - xx[0][i]) / 2;
      a = (xx[1][i] + xx[0][i]) / 2 - xx[2][i];
      ndnp->x[i+3] = (2 * a + b) / speed_intv;
    } else
      ndnp->x[i+3] = 0;
    ndnp->teval = tjd;
    ndnp->iephe = epheflag;
  }

  ndap = &swed.nddat[SEI_OSCU_APOG];
  Gmsm = GEOGCONST * (1 + 1 / EARTH_MOON_MRAT) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;

  for (i = istart; i <= 2; i++) {

    rxy =  sqrt(xx[i][0] * xx[i][0] + xx[i][1] * xx[i][1]);
    cosnode = xx[i][0] / rxy;
    sinnode = xx[i][1] / rxy;

    swi_cross_prod(xpos[i], xpos[i]+3, xnorm);
    rxy =  xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
    c2 = (rxy + xnorm[2] * xnorm[2]);
    rxyz = sqrt(c2);
    rxy = sqrt(rxy);
    sinincl = rxy / rxyz;
    cosincl = sqrt(1 - sinincl * sinincl);

    cosu = xpos[i][0] * cosnode + xpos[i][1] * sinnode;
    sinu = xpos[i][2] / sinincl;
    uu = atan2(sinu, cosu);

    rxyz = sqrt(square_sum(xpos[i]));
    v2 = square_sum((xpos[i]+3));
    sema = 1 / (2 / rxyz - v2 / Gmsm);

    pp = c2 / Gmsm;
    ecce = sqrt(1 - pp / sema);

    cosE = 1 / ecce * (1 - rxyz / sema);
    sinE = 1 / ecce / sqrt(sema * Gmsm) * dot_prod(xpos[i], (xpos[i]+3));

    ny = 2 * atan(sqrt((1+ecce)/(1-ecce)) * sinE / (1 + cosE));

    xxa[i][0] = swi_mod2PI(uu - ny + PI);
    xxa[i][1] = 0;
    xxa[i][2] = sema * (1 + ecce);

    swi_polcart(xxa[i], xxa[i]);
    swi_coortrf2(xxa[i], xxa[i], -sinincl, cosincl);
    swi_cartpol(xxa[i], xxa[i]);

    xxa[i][0] += atan2(sinnode, cosnode);
    swi_polcart(xxa[i], xxa[i]);

    ny = swi_mod2PI(ny - uu);

    cosE = cos(2 * atan(tan(ny / 2) / sqrt((1+ecce) / (1-ecce))));

    r[0] = sema * (1 - ecce * cosE);

    r[1] = sqrt(square_sum(xx[i]));

    for (j = 0; j <= 2; j++)
      xx[i][j] *= r[0] / r[1];
  }

  for (i = 0; i <= 2; i++) {

    ndap->x[i] = xxa[2][i];
    if (iflag & SEFLG_SPEED) {
      ndap->x[i+3] = (xxa[1][i] - xxa[0][i]) / speed_intv / 2;
    } else {
      ndap->x[i+3] = 0;
    }
    ndap->teval = tjd;
    ndap->iephe = epheflag;

    ndnp->x[i] = xx[2][i];
    if (iflag & SEFLG_SPEED) {
      ndnp->x[i+3] = (xx[1][i] - xx[0][i]) / speed_intv / 2;
    } else {
      ndnp->x[i+3] = 0;
    }
  }

  for (j = 0; j <= 1; j++) {
    double x[6];
    if (j == 0) {
      ndp = &swed.nddat[SEI_TRUE_NODE];
    } else {
      ndp = &swed.nddat[SEI_OSCU_APOG];
    }
    memset((void *) ndp->xreturn, 0, 24 * sizeof(double));

    for (i = 0; i <= 5; i++)
      ndp->xreturn[6+i] = ndp->x[i];

    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);

    swi_coortrf2(ndp->xreturn+6, ndp->xreturn+18, -oe->seps, oe->ceps);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(ndp->xreturn+9, ndp->xreturn+21, -oe->seps, oe->ceps);
#ifdef SID_TNODE_FROM_ECL_T0

    if (iflag & SEFLG_SIDEREAL) {

      swi_precess(ndp->xreturn+18, sip->t0, iflag, J_TO_J2000);
      if (iflag & SEFLG_SPEED)
	swi_precess_speed(ndp->xreturn+21, sip->t0, iflag, J_TO_J2000);
      if (!(iflag & SEFLG_J2000)) {

	swi_precess(ndp->xreturn+18, tjd, iflag, J2000_TO_J);
	if (iflag & SEFLG_SPEED)
	  swi_precess_speed(ndp->xreturn+21, tjd, iflag, J2000_TO_J);
      }
    }
#endif
    if (!(iflag & SEFLG_NONUT)) {
      swi_coortrf2(ndp->xreturn+18, ndp->xreturn+18, -swed.nut.snut, swed.nut.cnut);
      if (iflag & SEFLG_SPEED)
	swi_coortrf2(ndp->xreturn+21, ndp->xreturn+21, -swed.nut.snut, swed.nut.cnut);
    }

    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    ndp->xflgs = iflag;
    ndp->iephe = iflag & SEFLG_EPHMASK;
#ifdef SID_TNODE_FROM_ECL_T0

#else
    if (iflag & SEFLG_SIDEREAL) {

      if ((swed.sidd.sid_mode & SE_SIDBIT_ECL_T0)
        || (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE)) {
	for (i = 0; i <= 5; i++)
	  x[i] = ndp->xreturn[18+i];

	if (!(iflag & SEFLG_NONUT))
	  swi_nutate(x, iflag, TRUE);

	swi_precess(x, tjd, iflag, J_TO_J2000);
	if (iflag & SEFLG_SPEED)
	  swi_precess_speed(x, tjd, iflag, J_TO_J2000);
        if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
	  swi_trop_ra2sid_lon(x, ndp->xreturn+6, ndp->xreturn+18, iflag);

        } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
          swi_trop_ra2sid_lon_sosy(x, ndp->xreturn+6, iflag);
	}

	swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
        swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);

      } else {
	swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
	if (swi_get_ayanamsa_with_speed(ndp->teval, iflag, daya, serr) == ERR)
	  return ERR;
	ndp->xreturn[0] -= daya[0] * DEGTORAD;
	ndp->xreturn[3] -= daya[1] * DEGTORAD;
	swi_polcart_sp(ndp->xreturn, ndp->xreturn+6);
      }
    } else if (iflag & SEFLG_J2000) {

      for (i = 0; i <= 5; i++)
        x[i] = ndp->xreturn[18+i];

      swi_precess(x, tjd, iflag, J_TO_J2000);
      if (iflag & SEFLG_SPEED)
        swi_precess_speed(x, tjd, iflag, J_TO_J2000);
      for (i = 0; i <= 5; i++)
        ndp->xreturn[18+i] = x[i];
      swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
      swi_coortrf2(ndp->xreturn+18, ndp->xreturn+6, swed.oec2000.seps, swed.oec2000.ceps);
      if (iflag & SEFLG_SPEED)
        swi_coortrf2(ndp->xreturn+21, ndp->xreturn+9, swed.oec2000.seps, swed.oec2000.ceps);
      swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
    }
#endif

      for (i = 0; i < 2; i++) {
        ndp->xreturn[i] *= RADTODEG;
        ndp->xreturn[i+3] *= RADTODEG;
        ndp->xreturn[i+12] *= RADTODEG;
        ndp->xreturn[i+15] *= RADTODEG;
      }
      ndp->xreturn[0] = swe_degnorm(ndp->xreturn[0]);
      ndp->xreturn[12] = swe_degnorm(ndp->xreturn[12]);

  }
  return OK;
}

static int intp_apsides(double tjd, int ipl, int32 iflag, char *serr)
{
  int i;
  int32 flg1, flg2;
  struct plan_data *ndp;
  struct epsilon *oe;
  struct nut *nut;
  double daya[2];
  double speed_intv = 0.1;
  double t, dt;
  double xpos[3][6], xx[6], x[6];
  int32 speedf1, speedf2;
  oe = &swed.oec;
  nut = &swed.nut;
  ndp = &swed.nddat[ipl];

  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = ndp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  speedf1 = ndp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == ndp->teval
	&& tjd != 0
	&& flg1 == flg2
	&& (!speedf2 || speedf1)) {
    ndp->xflgs = iflag;
    ndp->iephe = iflag & SEFLG_MOSEPH;
    return OK;
  }

  for (t = tjd - speed_intv, i = 0; i < 3; t += speed_intv, i++) {
    if (! (iflag & SEFLG_SPEED) && i != 1) continue;
    swi_intp_apsides(t, xpos[i], ipl);
  }

  for (i = 0; i < 3; i++) {
    xx[i] = xpos[1][i];
    xx[i+3] = 0;
  }
  if (iflag & SEFLG_SPEED) {
    xx[3] = swe_difrad2n(xpos[2][0], xpos[0][0]) / speed_intv / 2.0;
    xx[4] = (xpos[2][1] - xpos[0][1]) / speed_intv / 2.0;
    xx[5] = (xpos[2][2] - xpos[0][2]) / speed_intv / 2.0;
  }
  memset((void *) ndp->xreturn, 0, 24 * sizeof(double));

  swi_polcart_sp(xx, xx);

  if (!(iflag & SEFLG_TRUEPOS)) {
    dt = sqrt(square_sum(xx)) * AUNIT / CLIGHT / 86400.0;
    for (i = 1; i < 3; i++)
      xx[i] -= dt * xx[i+3];
  }
  for (i = 0; i <= 5; i++)
    ndp->xreturn[i+6] = xx[i];

  swi_coortrf2(ndp->xreturn+6, ndp->xreturn+18, -oe->seps, oe->ceps);
  if (iflag & SEFLG_SPEED)
    swi_coortrf2(ndp->xreturn+9, ndp->xreturn+21, -oe->seps, oe->ceps);
  ndp->teval = tjd;
  ndp->xflgs = iflag;
  ndp->iephe = iflag & SEFLG_EPHMASK;
  if (iflag & SEFLG_SIDEREAL) {

    if ((swed.sidd.sid_mode & SE_SIDBIT_ECL_T0)
	|| (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE)) {
      for (i = 0; i <= 5; i++)
	x[i] = ndp->xreturn[18+i];

      swi_precess(x, tjd, iflag, J_TO_J2000);
      if (iflag & SEFLG_SPEED)
	swi_precess_speed(x, tjd, iflag, J_TO_J2000);
      if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
	swi_trop_ra2sid_lon(x, ndp->xreturn+6, ndp->xreturn+18, iflag);

      } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
	swi_trop_ra2sid_lon_sosy(x, ndp->xreturn+6, iflag);
      }

      swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
      swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    } else {

      swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
      if (swi_get_ayanamsa_with_speed(ndp->teval, iflag, daya, serr) == ERR)
        return ERR;
      ndp->xreturn[0] -= daya[0] * DEGTORAD;
      ndp->xreturn[3] -= daya[1] * DEGTORAD;
      swi_polcart_sp(ndp->xreturn, ndp->xreturn+6);
      swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    }
  } else if (iflag & SEFLG_J2000) {

    for (i = 0; i <= 5; i++)
      x[i] = ndp->xreturn[18+i];

    swi_precess(x, tjd, iflag, J_TO_J2000);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(x, tjd, iflag, J_TO_J2000);
    for (i = 0; i <= 5; i++)
      ndp->xreturn[18+i] = x[i];
    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    swi_coortrf2(ndp->xreturn+18, ndp->xreturn+6, swed.oec2000.seps, swed.oec2000.ceps);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(ndp->xreturn+21, ndp->xreturn+9, swed.oec2000.seps, swed.oec2000.ceps);
    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
  } else {

    if (!(iflag & SEFLG_NONUT)) {
      swi_nutate(ndp->xreturn+18, iflag, FALSE);
    }

    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);

    swi_coortrf2(ndp->xreturn+18, ndp->xreturn+6, oe->seps, oe->ceps);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(ndp->xreturn+21, ndp->xreturn+9, oe->seps, oe->ceps);
    if (!(iflag & SEFLG_NONUT)) {
      swi_coortrf2(ndp->xreturn+6, ndp->xreturn+6, nut->snut, nut->cnut);
      if (iflag & SEFLG_SPEED)
	swi_coortrf2(ndp->xreturn+9, ndp->xreturn+9, nut->snut, nut->cnut);
    }

    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
  }

  for (i = 0; i < 2; i++) {
    ndp->xreturn[i] *= RADTODEG;
    ndp->xreturn[i+3] *= RADTODEG;
    ndp->xreturn[i+12] *= RADTODEG;
    ndp->xreturn[i+15] *= RADTODEG;
  }
  ndp->xreturn[0] = swe_degnorm(ndp->xreturn[0]);
  ndp->xreturn[12] = swe_degnorm(ndp->xreturn[12]);

  return OK;
}

int swi_plan_for_osc_elem(int32 iflag, double tjd, double *xx)
{
  int i;
  double x[6];
  struct nut nuttmp;
  struct nut *nutp = &nuttmp;
  struct epsilon *oe = &swed.oec;
  struct epsilon oectmp;

  if (!(iflag & SEFLG_ICRS) && swi_get_denum(SEI_SUN, iflag) >= 403) {
    swi_bias(xx, tjd, iflag, FALSE);
  }

#ifdef SID_TNODE_FROM_ECL_T0
  struct sid_data *sip = &swed.sidd;

  if (iflag & SEFLG_SIDEREAL) {
    tjd = sip->t0;
    swi_precess(xx, tjd, iflag, J2000_TO_J);
    swi_precess(xx+3, tjd, iflag, J2000_TO_J);
    calc_epsilon(tjd, iflag, &oectmp);
    oe = &oectmp;
  } else if (!(iflag & SEFLG_J2000)) {
#endif
    swi_precess(xx, tjd, iflag, J2000_TO_J);
    swi_precess(xx+3, tjd, iflag, J2000_TO_J);

    if (tjd == swed.oec.teps) {
      oe = &swed.oec;
    } else if (tjd == J2000) {
      oe = &swed.oec2000;
    } else {
      calc_epsilon(tjd, iflag, &oectmp);
      oe = &oectmp;
    }
#ifdef SID_TNODE_FROM_ECL_T0
  } else {
    oe = &swed.oec2000;
  }
#endif

  if (!(iflag & SEFLG_NONUT)) {
    if (tjd == swed.nut.tnut) {
      nutp = &swed.nut;
    } else if (tjd == J2000) {
      nutp = &swed.nut2000;
    } else if (tjd == swed.nutv.tnut) {
      nutp = &swed.nutv;
    } else {
      nutp = &nuttmp;
      swi_nutation(tjd, iflag, nutp->nutlo);
      nutp->tnut = tjd;
      nutp->snut = sin(nutp->nutlo[1]);
      nutp->cnut = cos(nutp->nutlo[1]);
      nut_matrix(nutp, oe);
    }
    for (i = 0; i <= 2; i++) {
      x[i] = xx[0] * nutp->matrix[0][i] +
	     xx[1] * nutp->matrix[1][i] +
	     xx[2] * nutp->matrix[2][i];
    }

    for (i = 0; i <= 2; i++) {
      x[i+3] = xx[3] * nutp->matrix[0][i] +
	       xx[4] * nutp->matrix[1][i] +
	       xx[5] * nutp->matrix[2][i];
    }
    for (i = 0; i <= 5; i++)
      xx[i] = x[i];
  }

  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
#ifdef SID_TNODE_FROM_ECL_T0
  if (iflag & SEFLG_SIDEREAL) {

    swi_cartpol_sp(xx, xx);
    xx[0] -= sip->ayan_t0;
    swi_polcart_sp(xx, xx);
  } else
#endif
  if (!(iflag & SEFLG_NONUT)) {
    swi_coortrf2(xx, xx, nutp->snut, nutp->cnut);
    swi_coortrf2(xx+3, xx+3, nutp->snut, nutp->cnut);
  }
  return(OK);
}

static const struct meff_ele eff_arr[] = {

  {1.000, 1.000000},
  {0.990, 0.999979},
  {0.980, 0.999940},
  {0.970, 0.999881},
  {0.960, 0.999811},
  {0.950, 0.999724},
  {0.940, 0.999622},
  {0.930, 0.999497},
  {0.920, 0.999354},
  {0.910, 0.999192},
  {0.900, 0.999000},
  {0.890, 0.998786},
  {0.880, 0.998535},
  {0.870, 0.998242},
  {0.860, 0.997919},
  {0.850, 0.997571},
  {0.840, 0.997198},
  {0.830, 0.996792},
  {0.820, 0.996316},
  {0.810, 0.995791},
  {0.800, 0.995226},
  {0.790, 0.994625},
  {0.780, 0.993991},
  {0.770, 0.993326},
  {0.760, 0.992598},
  {0.750, 0.991770},
  {0.740, 0.990873},
  {0.730, 0.989919},
  {0.720, 0.988912},
  {0.710, 0.987856},
  {0.700, 0.986755},
  {0.690, 0.985610},
  {0.680, 0.984398},
  {0.670, 0.982986},
  {0.660, 0.981437},
  {0.650, 0.979779},
  {0.640, 0.978024},
  {0.630, 0.976182},
  {0.620, 0.974256},
  {0.610, 0.972253},
  {0.600, 0.970174},
  {0.590, 0.968024},
  {0.580, 0.965594},
  {0.570, 0.962797},
  {0.560, 0.959758},
  {0.550, 0.956515},
  {0.540, 0.953088},
  {0.530, 0.949495},
  {0.520, 0.945741},
  {0.510, 0.941838},
  {0.500, 0.937790},
  {0.490, 0.933563},
  {0.480, 0.928668},
  {0.470, 0.923288},
  {0.460, 0.917527},
  {0.450, 0.911432},
  {0.440, 0.905035},
  {0.430, 0.898353},
  {0.420, 0.891022},
  {0.410, 0.882940},
  {0.400, 0.874312},
  {0.390, 0.865206},
  {0.380, 0.855423},
  {0.370, 0.844619},
  {0.360, 0.833074},
  {0.350, 0.820876},
  {0.340, 0.808031},
  {0.330, 0.793962},
  {0.320, 0.778931},
  {0.310, 0.763021},
  {0.300, 0.745815},
  {0.290, 0.727557},
  {0.280, 0.708234},
  {0.270, 0.687583},
  {0.260, 0.665741},
  {0.250, 0.642597},
  {0.240, 0.618252},
  {0.230, 0.592586},
  {0.220, 0.565747},
  {0.210, 0.537697},
  {0.200, 0.508554},
  {0.190, 0.478420},
  {0.180, 0.447322},
  {0.170, 0.415454},
  {0.160, 0.382892},
  {0.150, 0.349955},
  {0.140, 0.316691},
  {0.130, 0.283565},
  {0.120, 0.250431},
  {0.110, 0.218327},
  {0.100, 0.186794},
  {0.090, 0.156287},
  {0.080, 0.128421},
  {0.070, 0.102237},
  {0.060, 0.077393},
  {0.050, 0.054833},
  {0.040, 0.036361},
  {0.030, 0.020953},
  {0.020, 0.009645},
  {0.010, 0.002767},
  {0.000, 0.000000}
};
static double meff(double r)
{
  double f, m;
  int i;
  if (r <= 0) {
    return 0.0;
  } else if (r >= 1) {
    return 1.0;
  }
  for (i = 0; eff_arr[i].r > r; i++)
    ;
  f = (r - eff_arr[i-1].r) / (eff_arr[i].r - eff_arr[i-1].r);
  m = eff_arr[i-1].m + f * (eff_arr[i].m - eff_arr[i-1].m);
  return m;
}

static void denormalize_positions(double *x0, double *x1, double *x2)
{
  int i;

  for (i = 0; i <= 12; i += 12) {
    if (x1[i] - x0[i] < -180)
      x0[i] -= 360;
    if (x1[i] - x0[i] > 180)
      x0[i] += 360;
    if (x1[i] - x2[i] < -180)
      x2[i] -= 360;
    if (x1[i] - x2[i] > 180)
      x2[i] += 360;
  }
}

static void calc_speed(double *x0, double *x1, double *x2, double dt)
{
  int i, j, k;
  double a, b;
  for (j = 0; j <= 18; j += 6) {
    for (i = 0; i < 3; i++) {
      k = j + i;
      b = (x2[k] - x0[k]) / 2;
      a = (x2[k] + x0[k]) / 2 - x1[k];
      x1[k+3] = (2 * a + b) / dt;
    }
  }
}

void swi_check_ecliptic(double tjd, int32 iflag)
{
  if (swed.oec2000.teps != J2000) {
    calc_epsilon(J2000, iflag, &swed.oec2000);
  }
  if (tjd == J2000) {
    swed.oec.teps = swed.oec2000.teps;
    swed.oec.eps = swed.oec2000.eps;
    swed.oec.seps = swed.oec2000.seps;
    swed.oec.ceps = swed.oec2000.ceps;
    return;
  }
  if (swed.oec.teps != tjd || tjd == 0) {
    calc_epsilon(tjd, iflag, &swed.oec);
  }
}

void swi_check_nutation(double tjd, int32 iflag)
{
  int32 speedf1, speedf2;
  static TLS int32 nutflag = 0;
  double t;
  speedf1 = nutflag & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (!(iflag & SEFLG_NONUT)
	&& (tjd != swed.nut.tnut || tjd == 0
	|| (!speedf1 && speedf2))) {
    swi_nutation(tjd, iflag, swed.nut.nutlo);
    swed.nut.tnut = tjd;
    swed.nut.snut = sin(swed.nut.nutlo[1]);
    swed.nut.cnut = cos(swed.nut.nutlo[1]);
    nutflag = iflag;
    nut_matrix(&swed.nut, &swed.oec);
    if (iflag & SEFLG_SPEED) {

      t = tjd - NUT_SPEED_INTV;
      swi_nutation(t, iflag, swed.nutv.nutlo);
      swed.nutv.tnut = t;
      swed.nutv.snut = sin(swed.nutv.nutlo[1]);
      swed.nutv.cnut = cos(swed.nutv.nutlo[1]);
      nut_matrix(&swed.nutv, &swed.oec);
    }
  }
}

static int32 plaus_iflag(int32 iflag, int32 ipl, double tjd, char *serr)
{
  int32 epheflag = 0;
  int jplhor_model = swed.astro_models[SE_MODEL_JPLHOR_MODE];
  int jplhora_model = swed.astro_models[SE_MODEL_JPLHORA_MODE];
  if (jplhor_model == 0) jplhor_model = SEMOD_JPLHOR_DEFAULT;
  if (jplhora_model == 0) jplhora_model = SEMOD_JPLHORA_DEFAULT;

  if (iflag & SEFLG_JPLHOR)
    iflag &= ~SEFLG_JPLHOR_APPROX;

  if (iflag & SEFLG_TOPOCTR) {
    iflag = iflag & ~(SEFLG_HELCTR | SEFLG_BARYCTR);
  }

  if (iflag & SEFLG_BARYCTR)
    iflag = iflag & ~(SEFLG_HELCTR);
  if (iflag & SEFLG_HELCTR)
    iflag = iflag & ~(SEFLG_BARYCTR);

  if (iflag & (SEFLG_HELCTR|SEFLG_BARYCTR))
    iflag |= SEFLG_NOABERR | SEFLG_NOGDEFL;

  if (iflag & SEFLG_J2000)
    iflag |= SEFLG_NONUT;

  if (iflag & SEFLG_SIDEREAL) {
    iflag |= SEFLG_NONUT;
    iflag = iflag & ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);
  }

  if (iflag & SEFLG_TRUEPOS)
    iflag |= (SEFLG_NOGDEFL | SEFLG_NOABERR);
  if (iflag & SEFLG_MOSEPH)
    epheflag = SEFLG_MOSEPH;
  if (iflag & SEFLG_SWIEPH)
    epheflag = SEFLG_SWIEPH;
  if (iflag & SEFLG_JPLEPH)
    epheflag = SEFLG_JPLEPH;
  if (epheflag == 0)
    epheflag = SEFLG_DEFAULTEPH;
  iflag = (iflag & ~SEFLG_EPHMASK) | epheflag;

  if (!(epheflag & SEFLG_JPLEPH))
    iflag = iflag & ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);

  if (ipl == SE_OSCU_APOG || ipl == SE_TRUE_NODE
      || ipl == SE_MEAN_APOG || ipl == SE_MEAN_NODE
      || ipl == SE_INTP_APOG || ipl == SE_INTP_PERG)
    iflag = iflag & ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);
  if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX)
    iflag = iflag & ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);

  if (iflag & SEFLG_JPLHOR) {
    if (swed.eop_dpsi_loaded <= 0) {
      if (serr != NULL) {
	switch (swed.eop_dpsi_loaded) {
	  case 0:
	    strcpy(serr, "you did not call swe_set_jpl_file(); default to SEFLG_JPLHOR_APPROX");
	    break;
	  case -1:
	    strcpy(serr, "file eop_1962_today.txt not found; default to SEFLG_JPLHOR_APPROX");
	    break;
	  case -2:
	    strcpy(serr, "file eop_1962_today.txt corrupt; default to SEFLG_JPLHOR_APPROX");
	    break;
	  case -3:
	    strcpy(serr, "file eop_finals.txt corrupt; default to SEFLG_JPLHOR_APPROX");
	    break;
	}
      }
      iflag &= ~SEFLG_JPLHOR;
      iflag |= SEFLG_JPLHOR_APPROX;
    }
  }
  if (iflag & SEFLG_JPLHOR)
    iflag |= SEFLG_ICRS;
  if ((iflag & SEFLG_JPLHOR_APPROX) && jplhora_model == SEMOD_JPLHORA_2)
    iflag |= SEFLG_ICRS;
  return iflag;
}

static int32 fixstar_format_search_name(char *star, char *sstar, char *serr)
{
  char *sp;
  size_t cmplen;
  strncpy(sstar, star, SWI_STAR_LENGTH);
  sstar[SWI_STAR_LENGTH] = '\0';

  while ((sp = strchr(sstar, ' ')) != NULL)
    swi_strcpy(sp, sp+1);

  for (sp = sstar; *sp != '\0' && *sp != ','; sp++)
    *sp = tolower((int) *sp);
  cmplen = strlen(sstar);
  if (cmplen == 0) {
    if (serr != NULL)
      sprintf(serr, "swe_fixstar(): star name empty");
    return ERR;
  }
  return OK;
}

static int32 save_star_in_struct(int nrecs, struct fixed_star *fstp, char *serr)
{
  int sizestru = sizeof(struct fixed_star);
  struct fixed_star *ftarget;
  char *serr_alloc = "error in function load_all_fixed_stars(): could not resize fixed stars array";
  if ((swed.fixed_stars = (struct fixed_star *) realloc(swed.fixed_stars, nrecs * sizestru)) == NULL) {
    if (serr != NULL) strcpy(serr, serr_alloc);
    return ERR;
  }
  ftarget = swed.fixed_stars + (nrecs - 1);
  memcpy((void *) ftarget, (void *) fstp, sizestru);
  return OK;
}

static int CMP_CALL_CONV fixedstar_name_compare(const void *star1, const void *star2)
{
  const struct fixed_star *fs1 = (const struct fixed_star *) star1;
  const struct fixed_star *fs2 = (const struct fixed_star *) star2;
  return strcmp(fs1->skey, fs2->skey);
}

static int CMP_CALL_CONV fstar_node_compare(const void *node1, const void *node2)
{
  const struct fixed_star *n1 = (const struct fixed_star *) node1;
  const struct fixed_star *n2 = (const struct fixed_star *) node2;
  return strcmp(n1->skey, n2->skey);
}

int32 fixstar_cut_string(char *srecord, char *star, struct fixed_star *stardata, char *serr)
{
  int i;
  char s[AS_MAXCH + 20];
  char *sde_d;
  char *cpos[20];
  double epoch, radv, parall, mag;
  double ra_s, ra_pm, de_pm, ra, de;
  double ra_h, ra_m, de_d, de_m, de_s;
  strcpy(s, srecord);
  i = swi_cutstr(s, ",", cpos, 20);

  swi_right_trim(cpos[0]);
  swi_right_trim(cpos[1]);
  if (i < 14) {
    if (serr != NULL) {
      if (i >= 2) {
	sprintf(serr, "data of star '%s,%s' incomplete", cpos[0], cpos[1]);
      } else {
        if (strlen(s) > 200) s[200] = '\0';
	sprintf(serr, "invalid line in fixed stars file: '%s'", s);
      }
    }
    return ERR;
  }
  if (strlen(cpos[0]) > SWI_STAR_LENGTH)
    cpos[0][SWI_STAR_LENGTH] = '\0';
  if (strlen(cpos[1]) > SWI_STAR_LENGTH-1)
    cpos[1][SWI_STAR_LENGTH-1] = '\0';
  if (star != NULL) {
    strcpy(star, cpos[0]);
    if (strlen(cpos[0]) + strlen(cpos[1]) + 1 < SWI_STAR_LENGTH - 1)
      sprintf(star + strlen(star), ",%s", cpos[1]);
  }
  strcpy(stardata->starname, cpos[0]);
  strcpy(stardata->starbayer, cpos[1]);

  epoch = atof(cpos[2]);
  ra_h = atof(cpos[3]);
  ra_m = atof(cpos[4]);
  ra_s = atof(cpos[5]);
  de_d = atof(cpos[6]);
  sde_d = cpos[6];
  de_m = atof(cpos[7]);
  de_s = atof(cpos[8]);
  ra_pm = atof(cpos[9]);
  de_pm = atof(cpos[10]);
  radv = atof(cpos[11]);
  parall = atof(cpos[12]);
  if (parall < 0) parall = -parall;
  mag = atof(cpos[13]);

  ra = (ra_s / 3600.0 + ra_m / 60.0 + ra_h) * 15.0;
  if (strchr(sde_d, '-') == NULL) {
    de = de_s / 3600.0 + de_m / 60.0 + de_d;
  } else {
    de = -de_s / 3600.0 - de_m / 60.0 + de_d;
  }

  if (swed.is_old_starfile == TRUE) {
    ra_pm = ra_pm * 15 / 3600.0;
    de_pm = de_pm / 3600.0;
  } else {
    ra_pm = ra_pm / 10.0 / 3600.0;
    de_pm = de_pm / 10.0 / 3600.0;
    parall /= 1000.0;
  }

  if (parall > 1) {
    parall = (1 / parall / 3600.0);
  } else {
    parall /= 3600;
  }

  radv *= KM_S_TO_AU_CTY;

  ra *= DEGTORAD;
  de *= DEGTORAD;
  ra_pm *= DEGTORAD;
  de_pm *= DEGTORAD;
  ra_pm /= cos(de);
  parall *= DEGTORAD;
  stardata->epoch = epoch;
  stardata->ra = ra;
  stardata->de = de;
  stardata->ramot = ra_pm;
  stardata->demot = de_pm;
  stardata->parall = parall;
  stardata->radvel = radv;
  stardata->mag = mag;
  return OK;
}

static int32 load_all_fixed_stars(char *serr)
{
  int32 retc = OK;
  int nstars = 0, nrecs = 0, nnamed = 0;
  char s[AS_MAXCH], *sp;
  char srecord[AS_MAXCH];
  struct fixed_star fstdata;
  char last_starbayer[SWI_STAR_LENGTH + 1];
  *last_starbayer = '\0';
  if (swed.n_fixstars_records > 0) {
    return -2;
  }
  if (swed.fixfp == NULL) {
    if ((swed.fixfp = swi_fopen(SEI_FILE_FIXSTAR, SE_STARFILE, swed.ephepath, serr)) == NULL) {
      swed.is_old_starfile = TRUE;
      if ((swed.fixfp = swi_fopen(SEI_FILE_FIXSTAR, SE_STARFILE_OLD, swed.ephepath, NULL)) == NULL) {
	swed.is_old_starfile = FALSE;

	return ERR;
      }
    }
  }
  rewind(swed.fixfp);
  swed.fixed_stars = NULL;
  while (fgets(s, AS_MAXCH, swed.fixfp) != NULL) {

    if (*s == '#') continue;
    if (*s == '\n') continue;
    if (*s == '\r') continue;
    if (*s == '\0') continue;
    strcpy(srecord, s);
    retc = fixstar_cut_string(srecord, NULL, &fstdata, serr);
    if (retc == ERR) return ERR;

    if (*fstdata.starname != '\0') {
      nrecs++;
      nnamed++;
      strcpy(fstdata.skey, fstdata.starname);

      while ((sp = strchr(fstdata.skey, ' ')) != NULL)
	swi_strcpy(sp, sp+1);

      for (sp = fstdata.skey; *sp != '\0'; sp++)
	*sp = tolower((int) *sp);
      if ((retc = save_star_in_struct(nrecs, &fstdata, serr)) == ERR) return ERR;
    }

    if (strcmp(fstdata.starbayer, last_starbayer) == 0)
      continue;
    nstars++;
    nrecs++;

    sprintf(fstdata.skey, ",%s", fstdata.starbayer);

    while ((sp = strchr(fstdata.skey, ' ')) != NULL)
      swi_strcpy(sp, sp+1);
    strcpy(last_starbayer, fstdata.starbayer);
    if ((retc = save_star_in_struct(nrecs, &fstdata, serr)) == ERR) return ERR;

  }
  swed.n_fixstars_real = nstars;
  swed.n_fixstars_named = nnamed;
  swed.n_fixstars_records = nrecs;

  (void) qsort ((void *) swed.fixed_stars, (size_t) nrecs, sizeof (struct fixed_star),
                    (int (CMP_CALL_CONV *)(const void *,const void *))(fixedstar_name_compare));
  return retc;
}

static int32 fixstar_calc_from_struct(struct fixed_star *stardata, double tjd, int32 iflag, char *star, double *xx, char *serr)
{
  int i;
  int32 retc = OK;
  double epoch, radv, parall;
  double ra_pm, de_pm, ra, de, t;
  double daya[2], rdist;
  double x[6], xxsv[6], xobs[6], xobs_dt[6], *xpo = NULL, *xpo_dt = NULL;
  static TLS double xearth[6], xearth_dt[6], xsun[6], xsun_dt[6];
  double dt = PLAN_SPEED_INTV * 0.1;
  int32 epheflag, iflgsave;
  struct epsilon *oe = &swed.oec2000;
  iflgsave = iflag;
  iflag |= SEFLG_SPEED;
  if (serr != NULL)
    *serr = '\0';
  iflag = plaus_iflag(iflag, -1, tjd, serr);
  epheflag = iflag & SEFLG_EPHMASK;
  if (swi_init_swed_if_start() == 1 && !(epheflag & SEFLG_MOSEPH) && serr != NULL) {
    strcpy(serr, "Please call swe_set_ephe_path() or swe_set_jplfile() before calling swe_fixstar() or swe_fixstar_ut()");
  }
  if (swed.last_epheflag != epheflag) {
    free_planets();

    if (swed.jpl_file_is_open) {
      swi_close_jpl_file();
      swed.jpl_file_is_open = FALSE;
    }
    for (i = 0; i < SEI_NEPHFILES; i ++) {
      if (swed.fidat[i].fptr != NULL)
	fclose(swed.fidat[i].fptr);
      memset((void *) &swed.fidat[i], 0, sizeof(struct file_data));
    }
    swed.last_epheflag = epheflag;
  }

  if (iflag & SEFLG_SIDEREAL && !swed.ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);

  swi_check_ecliptic(tjd, iflag);

  swi_check_nutation(tjd, iflag);
  sprintf(star, "%s,%s", stardata->starname, stardata->starbayer);
  epoch = stardata->epoch;
  ra_pm = stardata->ramot; de_pm = stardata->demot;
  radv = stardata->radvel; parall = stardata->parall;
  ra = stardata->ra; de = stardata->de;
  if (epoch == 1950) {
    t= (tjd - B1950);
  } else {
    t= (tjd - J2000);
  }
  x[0] = ra;
  x[1] = de;
  x[2] = 1;
  if (parall == 0) {
    rdist = 1000000000;
  } else {
    rdist = 1.0 / (parall * RADTODEG * 3600) * PARSEC_TO_AUNIT;

  }

  x[2] = rdist;
  x[3] = ra_pm / 36525.0;
  x[4] = de_pm / 36525.0;
  x[5] = radv / 36525.0;

  swi_polcart_sp(x, x);

  if (epoch == 1950) {
    swi_FK4_FK5(x, B1950);
    swi_precess(x, B1950, 0, J_TO_J2000);
    swi_precess(x+3, B1950, 0, J_TO_J2000);
  }

  if (epoch != 0) {
    swi_icrs2fk5(x, iflag, TRUE);

    if (swi_get_denum(SEI_SUN, iflag) >= 403) {
      swi_bias(x, J2000, SEFLG_SPEED, FALSE);
    }
  }

  if (!(iflag & SEFLG_BARYCTR) && (!(iflag & SEFLG_HELCTR) || !(iflag & SEFLG_MOSEPH))) {
    if ((retc =  main_planet_bary(tjd - dt, SEI_EARTH, epheflag, iflag, NO_SAVE, xearth_dt, xearth_dt, xsun_dt, NULL, serr)) != OK) {
      return ERR;
    }
    if ((retc =  main_planet_bary(tjd, SEI_EARTH, epheflag, iflag, DO_SAVE, xearth, xearth, xsun, NULL, serr)) != OK) {
      return ERR;
    }
  }

  if (iflag & SEFLG_TOPOCTR) {
    if (swi_get_observer(tjd - dt, iflag | SEFLG_NONUT, NO_SAVE, xobs_dt, serr) != OK)
      return ERR;
    if (swi_get_observer(tjd, iflag | SEFLG_NONUT, NO_SAVE, xobs, serr) != OK)
      return ERR;

    for (i = 0; i <= 5; i++) {
      xobs[i] = xobs[i] + xearth[i];
      xobs_dt[i] = xobs_dt[i] + xearth_dt[i];
    }
  } else if (!(iflag & SEFLG_BARYCTR) && (!(iflag & SEFLG_HELCTR) || !(iflag & SEFLG_MOSEPH))) {

    for (i = 0; i <= 5; i++) {
      xobs[i] = xearth[i];
      xobs_dt[i] = xearth_dt[i];
    }
  }

  if ((iflag & SEFLG_HELCTR) && (iflag & SEFLG_MOSEPH)) {
    xpo = NULL;
    xpo_dt = NULL;
  } else if (iflag & SEFLG_HELCTR) {
    xpo = xsun;
    xpo_dt = xsun_dt;
  } else if (iflag & SEFLG_BARYCTR) {
    xpo = NULL;
    xpo_dt = NULL;
  } else {
    xpo = xobs;
    xpo_dt = xobs_dt;
  }
  if (xpo == NULL) {
    for (i = 0; i <= 2; i++) {
      x[i] += t * x[i+3];
    }
  } else {
    for (i = 0; i <= 2; i++) {
      x[i] += t * x[i+3];
      x[i] -= xpo[i];
      x[i+3] -= xpo[i+3];
    }
  }

  if ((iflag & SEFLG_TRUEPOS) == 0 && (iflag & SEFLG_NOGDEFL) == 0) {
    swi_deflect_light(x, 0, iflag & SEFLG_SPEED);
  }

  if ((iflag & SEFLG_TRUEPOS) == 0 && (iflag & SEFLG_NOABERR) == 0)
    swi_aberr_light_ex(x, xpo, xpo_dt, dt, iflag & SEFLG_SPEED);

  if (!(iflag & SEFLG_ICRS) && (swi_get_denum(SEI_SUN, iflag) >= 403 || (iflag & SEFLG_BARYCTR))) {
    swi_bias(x, tjd, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = x[i];

  if ((iflag & SEFLG_J2000) == 0) {
    swi_precess(x, tjd, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(x, tjd, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else
    oe = &swed.oec2000;

  if (!(iflag & SEFLG_NONUT))
    swi_nutate(x, iflag, FALSE);
if ((0)) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  printf("%.17f %.17f %f\n", x[0]/r, x[1]/r, x[2]/r);
}

  if ((iflag & SEFLG_EQUATORIAL) == 0) {
    swi_coortrf2(x, x, oe->seps, oe->ceps);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(x+3, x+3, oe->seps, oe->ceps);
    if (!(iflag & SEFLG_NONUT)) {
      swi_coortrf2(x, x, swed.nut.snut, swed.nut.cnut);
      if (iflag & SEFLG_SPEED)
	swi_coortrf2(x+3, x+3, swed.nut.snut, swed.nut.cnut);
    }
  }

  if (iflag & SEFLG_SIDEREAL) {

    if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
      if (swi_trop_ra2sid_lon(xxsv, x, xxsv, iflag) != OK)
        return ERR;
      if (iflag & SEFLG_EQUATORIAL) {
        for (i = 0; i <= 5; i++)
          x[i] = xxsv[i];
      }

    } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
      if (swi_trop_ra2sid_lon_sosy(xxsv, x, iflag) != OK)
	return ERR;
      if (iflag & SEFLG_EQUATORIAL) {
        for (i = 0; i <= 5; i++)
          x[i] = xxsv[i];
      }

    } else {
      swi_cartpol_sp(x, x);

      if (swi_get_ayanamsa_with_speed(tjd, iflag, daya, serr) == ERR)
        return ERR;
      x[0] -= daya[0] * DEGTORAD;
      x[3] -= daya[1] * DEGTORAD;
      swi_polcart_sp(x, x);
    }
  }

  if ((iflag & SEFLG_XYZ) == 0)
    swi_cartpol_sp(x, x);

  if ((iflag & SEFLG_RADIANS) == 0 && (iflag & SEFLG_XYZ) == 0) {
    for (i = 0; i < 2; i++) {
      x[i] *= RADTODEG;
      x[i+3] *= RADTODEG;
    }
  }
  for (i = 0; i <= 5; i++)
    xx[i] = x[i];
  if (!(iflgsave & SEFLG_SPEED)) {
    for (i = 3; i <= 5; i++)
      xx[i] = 0;
  }

  if ((iflgsave & SEFLG_EPHMASK) == 0)
    iflag = iflag & ~SEFLG_DEFAULTEPH;
  iflag = iflag & ~SEFLG_SPEED;
  return iflag;
}

static int32 search_star_in_list(char *sstar, struct fixed_star *stardata, char *serr)
{
  int i, star_nr = 0, ndata = 0, len;
  char *sp;
  char searchkey[AS_MAXCH];
  AS_BOOL is_bayer = FALSE;
  struct fixed_star *stardatap;
  struct fixed_star *stardatabegp;
  if (*sstar == ',') {
    is_bayer = TRUE;
  } else if (isdigit((int) *sstar)) {
    star_nr = atoi(sstar);
  } else {
    if ((sp = strchr(sstar, ',')) != NULL) {
      swi_strcpy(sstar, sp);
      is_bayer = TRUE;
    }
  }
  if (star_nr > 0) {
    if (star_nr > swed.n_fixstars_real) {
      if (serr != NULL)
	sprintf(serr, "error, swe_fixstar(): sequential fixed star number %d is not available", star_nr);
      return ERR;
    }
    *stardata = swed.fixed_stars[star_nr - 1];

    return OK;

  } else if (!is_bayer && (sp = strchr(sstar, '%')) != NULL) {
    stardatabegp = &(swed.fixed_stars[swed.n_fixstars_real]);
    ndata = swed.n_fixstars_named;
    if (sp - sstar != strlen(sstar) - 1) {
      if (serr != NULL)
	sprintf(serr, "error, swe_fixstar(): invalid search string %s", sstar);
      return ERR;
    }
    strcpy(searchkey, sstar);
    len = (int) (strlen(sstar) - 1);
    searchkey[len] = '\0';
    for (i = 0; i < ndata; i++) {
      if (strncmp(stardatabegp[i].skey, sstar, len) == 0) {
        *stardata = stardatabegp[i];
	return OK;
      }
    }
    if (serr != NULL)
      sprintf(serr, "error, swe_fixstar(): star search string %s did not match", sstar);
    return ERR;

  } else {
    strcpy(searchkey, sstar);
    if (is_bayer) {

      stardatabegp = swed.fixed_stars;
      ndata = swed.n_fixstars_real;
    } else {
      stardatabegp = &(swed.fixed_stars[swed.n_fixstars_real]);
      ndata = swed.n_fixstars_named;
    }
    stardatap = (struct fixed_star *) bsearch((void *) searchkey,
	       (void *) stardatabegp, (size_t) ndata,
	       sizeof (struct fixed_star),
	       fstar_node_compare);
    if (stardatap == NULL) {
      if (serr != NULL)
	sprintf(serr, "error, swe_fixstar(): could not find star name %s", sstar);
      return ERR;
    }
    *stardata = *stardatap;

    return OK;
  }
}

static AS_BOOL get_builtin_star(char *star, char *sstar, char *srecord)
{

  if (strncmp(star, "spica", 5) == 0 || strncmp(star, "Spica", 5) == 0) {
    strcpy(srecord, "Spica,alVir,ICRS,13,25,11.57937,-11,09,40.7501,-42.35,-30.67,1,13.06,0.97,-10,3672");
    strcpy(sstar, "spica");
    return TRUE;

  } else if (strstr(star, ",zePsc") != NULL || strncmp(star, "revati", 6) == 0 || strncmp(star, "Revati", 6) == 0) {
    strcpy(srecord, "Revati,zePsc,ICRS,01,13,43.88735,+07,34,31.2745,145,-55.69,15,18.76,5.187,06,174");
    strcpy(sstar, "revati");
    return TRUE;

  } else if (strstr(star, ",deCnc") != NULL || strncmp(star, "pushya", 6) == 0 || strncmp(star, "Pushya", 6) == 0 ) {
    strcpy(srecord, "Pushya,deCnc,ICRS,08,44,41.09921,+18,09,15.5034,-17.67,-229.26,17.14,24.98,3.94,18,2027");
    strcpy(sstar, "pushya");
    return TRUE;

  } else if (strstr(star, ",deCnc") != NULL) {
    strcpy(srecord, "Pushya,deCnc,ICRS,08,44,41.09921,+18,09,15.5034,-17.67,-229.26,17.14,24.98,3.94,18,2027");
    strcpy(sstar, "pushya");
    return TRUE;

  } else if (strstr(star, ",laSco") != NULL || strncmp(star, "mula", 6) == 0 || strncmp(star, "Mula", 6) == 0) {
    strcpy(srecord, "Mula,laSco,ICRS,17,33,36.52012,-37,06,13.7648,-8.53,-30.8,-3,5.71,1.62,-37,11673");
    strcpy(sstar, "mula");
    return TRUE;

  } else if (strstr(star, ",SgrA*") != NULL) {
    strcpy(srecord, "Gal. Center,SgrA*,2000,17,45,40.03599,-29,00,28.1699,-2.755718425,-5.547,0.0,0.125,999.99,0,0");
    strcpy(sstar, ",SgrA*");
    return TRUE;

  } else if (strstr(star, ",GP1958") != NULL) {
    strcpy(srecord, "Gal. Pole IAU1958,GP1958,1950,12,49,0.0,27,24,0.0,0.0,0.0,0.0,0.0,0.0,0,0");
    strcpy(sstar, ",GP1958");
    return TRUE;

  } else if (strstr(star, ",GPol") != NULL) {
    strcpy(srecord, "Gal. Pole,GPol,ICRS,12,51,36.7151981,27,06,11.193172,0.0,0.0,0.0,0.0,0.0,0,0");
    strcpy(sstar, ",GPol");
    return TRUE;

  } else if (strstr(star, ",GPol") != NULL) {
    strcpy(srecord, "Gal. Pole,GPol,ICRS,12,51,36.7151981,27,06,11.193172,0.0,0.0,0.0,0.0,0.0,0,0");
    strcpy(sstar, ",GPol");
    return TRUE;
  }
  return FALSE;
}

int32 CALL_CONV swe_fixstar2(char *star, double tjd, int32 iflag,
  double *xx, char *serr)
{
  int i;
  AS_BOOL is_builtin_star = FALSE;
  char sstar[SWI_STAR_LENGTH + 1];

  static TLS char slast_starname[AS_MAXCH];
  static TLS struct fixed_star last_stardata;
  char srecord[AS_MAXCH + 20];
  int retc;
  struct fixed_star stardata;
  if (serr != NULL)
    *serr = '\0';
#ifdef TRACE
  swi_open_trace(serr);
  trace_swe_fixstar(1, star, tjd, iflag, xx, serr);
#endif
  load_all_fixed_stars(serr);
  retc = fixstar_format_search_name(star, sstar, serr);
  if (retc == ERR)
    goto return_err;

  if (swed.n_fixstars_records > 0 && strcmp(slast_starname, sstar) == 0) {

    stardata = last_stardata;
    goto found;
  }
  if (get_builtin_star(star, sstar, srecord)) {
    is_builtin_star = TRUE;
  }
  if (is_builtin_star) {
    retc = fixstar_cut_string(srecord, star, &stardata, serr);

    goto found;

  }
  retc = search_star_in_list(sstar, &stardata, serr);
  if (retc == ERR)
    goto return_err;

  found:

  last_stardata = stardata;
  strcpy(slast_starname, sstar);
  if ((retc = fixstar_calc_from_struct(&stardata, tjd, iflag, star, xx, serr)) == ERR)
    goto return_err;
#ifdef TRACE
  trace_swe_fixstar(2, star, tjd, iflag, xx, serr);
#endif
  return iflag;
  return_err:
  for (i = 0; i <= 5; i++)
    xx[i] = 0;
#ifdef TRACE
  trace_swe_fixstar(2, star, tjd, iflag, xx, serr);
#endif
  return retc;
}

int32 CALL_CONV swe_fixstar2_ut(char *star, double tjd_ut, int32 iflag,
  double *xx, char *serr)
{
  double deltat;
  int32 retflag;
  int32 epheflag = 0;
  iflag = plaus_iflag(iflag, -1, tjd_ut, serr);
  epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag == 0) {
    epheflag = SEFLG_SWIEPH;
    iflag |= SEFLG_SWIEPH;
  }
  deltat = swe_deltat_ex(tjd_ut, iflag, serr);

  retflag = swe_fixstar2(star, tjd_ut + deltat, iflag, xx, serr);
  if (retflag != ERR && (retflag & SEFLG_EPHMASK) != epheflag) {
    deltat = swe_deltat_ex(tjd_ut, retflag, NULL);
    retflag = swe_fixstar2(star, tjd_ut + deltat, iflag, xx, NULL);
  }
  return retflag;
}

int32 CALL_CONV swe_fixstar2_mag(char *star, double *mag, char *serr)
{
  char sstar[SWI_STAR_LENGTH + 1];

  static TLS char slast_starname[AS_MAXCH];
  static TLS struct fixed_star last_stardata;
  int retc;
  struct fixed_star stardata;
  if (serr != NULL)
    *serr = '\0';
  load_all_fixed_stars(serr);
  retc = fixstar_format_search_name(star, sstar, serr);
  if (retc == ERR)
    goto return_err;

  if (swed.n_fixstars_records > 0 && strcmp(slast_starname, sstar) == 0) {

    stardata = last_stardata;
    goto found;
  }
  retc = search_star_in_list(sstar, &stardata, serr);
  if (retc == ERR)
    goto return_err;

  found:
  last_stardata = stardata;
  strcpy(slast_starname, sstar);
  *mag = stardata.mag;
  sprintf(star, "%s,%s", stardata.starname, stardata.starbayer);
  return OK;
  return_err:
  *mag = 0;
  return retc;
}

char *CALL_CONV swe_get_planet_name(int ipl, char *s)
{
  int i;
  int32 retc;
  double xp[6];
#ifdef TRACE
  swi_open_trace(NULL);
  trace_swe_get_planet_name(1, ipl, s);
#endif
  swi_init_swed_if_start();

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  if (ipl != 0 && ipl == swed.i_saved_planet_name) {
    strcpy(s, swed.saved_planet_name);
    return s;
  }
  switch(ipl) {
    case SE_SUN:
      strcpy(s, SE_NAME_SUN);
      break;
    case SE_MOON:
      strcpy(s, SE_NAME_MOON);
      break;
    case SE_MERCURY:
      strcpy(s, SE_NAME_MERCURY);
      break;
    case SE_VENUS:
      strcpy(s, SE_NAME_VENUS);
      break;
    case SE_MARS:
      strcpy(s, SE_NAME_MARS);
      break;
    case SE_JUPITER:
      strcpy(s, SE_NAME_JUPITER);
      break;
    case SE_SATURN:
      strcpy(s, SE_NAME_SATURN);
      break;
    case SE_URANUS:
      strcpy(s, SE_NAME_URANUS);
      break;
    case SE_NEPTUNE:
      strcpy(s, SE_NAME_NEPTUNE);
      break;
    case SE_PLUTO:
      strcpy(s, SE_NAME_PLUTO);
      break;
    case SE_MEAN_NODE:
      strcpy(s, SE_NAME_MEAN_NODE);
      break;
    case SE_TRUE_NODE:
      strcpy(s, SE_NAME_TRUE_NODE);
      break;
    case SE_MEAN_APOG:
      strcpy(s, SE_NAME_MEAN_APOG);
      break;
    case SE_OSCU_APOG:
      strcpy(s, SE_NAME_OSCU_APOG);
      break;
    case SE_INTP_APOG:
      strcpy(s, SE_NAME_INTP_APOG);
      break;
    case SE_INTP_PERG:
      strcpy(s, SE_NAME_INTP_PERG);
      break;
    case SE_EARTH:
      strcpy(s, SE_NAME_EARTH);
      break;
    case SE_CHIRON:
    case SE_AST_OFFSET + MPC_CHIRON:
      strcpy(s, SE_NAME_CHIRON);
      break;
    case SE_PHOLUS:
    case SE_AST_OFFSET + MPC_PHOLUS:
      strcpy(s, SE_NAME_PHOLUS);
      break;
    case SE_CERES:
    case SE_AST_OFFSET + MPC_CERES:
      strcpy(s, SE_NAME_CERES);
      break;
    case SE_PALLAS:
    case SE_AST_OFFSET + MPC_PALLAS:
      strcpy(s, SE_NAME_PALLAS);
      break;
    case SE_JUNO:
    case SE_AST_OFFSET + MPC_JUNO:
      strcpy(s, SE_NAME_JUNO);
      break;
    case SE_VESTA:
    case SE_AST_OFFSET + MPC_VESTA:
      strcpy(s, SE_NAME_VESTA);
      break;
    default:

      if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX) {
        swi_get_fict_name(ipl - SE_FICT_OFFSET, s);
        break;
      }

      if (ipl > SE_PLMOON_OFFSET || ipl > SE_AST_OFFSET) {

	if (ipl == swed.fidat[SEI_FILE_ANY_AST].ipl[0]) {
	  strcpy(s, swed.fidat[SEI_FILE_ANY_AST].astnam);

	} else {
	  retc = sweph(J2000, ipl, SEI_FILE_ANY_AST, 0, NULL, NO_SAVE, xp, NULL);
	  if (retc != ERR && retc != NOT_AVAILABLE) {
	    strcpy(s, swed.fidat[SEI_FILE_ANY_AST].astnam);
	  } else {
	    if (ipl > SE_AST_OFFSET) {
	      sprintf(s, "%d: not found (asteroid)", ipl - SE_AST_OFFSET);
	    } else {
	      sprintf(s, "%d: not found (planetary moon)", ipl);
	    }
	  }
	}

        if (ipl > SE_AST_OFFSET && (s[0] == '?' || isdigit((int) s[1]))) {
          int ipli = (int) (ipl - SE_AST_OFFSET), iplf = 0;
          FILE *fp;
          char si[AS_MAXCH], *sp, *sp2;
          if ((fp = swi_fopen(-1, SE_ASTNAMFILE, swed.ephepath, NULL)) != NULL) {
            while(ipli != iplf && (sp = fgets(si, AS_MAXCH, fp)) != NULL) {
              while (*sp == ' ' || *sp == '\t'
                     || *sp == '(' || *sp == '[' || *sp == '{')
                sp++;
              if (*sp == '#' || *sp == '\r' || *sp == '\n' || *sp == '\0')
                continue;

              iplf = atoi(sp);
              if (ipli != iplf)
                continue;

              sp = strpbrk(sp, " \t");
              if (sp == NULL)
                continue;
              while (*sp == ' ' || *sp == '\t')
                sp++;
              sp2 = strpbrk(sp, "#\r\n");
              if (sp2 != NULL)
                *sp2 = '\0';
              if (*sp == '\0')
                continue;
              swi_right_trim(sp);
              strcpy(s, sp);
            }
            fclose(fp);
          }
        }
      } else  {
	i = ipl;
	sprintf(s, "%d", i);
      }
      break;
  }
#ifdef TRACE
  swi_open_trace(NULL);
  trace_swe_get_planet_name(2, ipl, s);
#endif
  if (strlen(s) < 80) {
    swed.i_saved_planet_name = ipl;
    strcpy(swed.saved_planet_name, s);
  }
  return s;
}

const char *CALL_CONV swe_get_ayanamsa_name(int32 isidmode)
{
  isidmode %= SE_SIDBITS;
  if (isidmode < SE_NSIDM_PREDEF)
    return ayanamsa_name[isidmode];
  return NULL;
}

#ifdef TRACE
static void trace_swe_calc(int swtch, double tjd, int ipl, int32 iflag, double *xx, char *serr)
{
  if (swi_trace_count >= TRACE_COUNT_MAX)
    return;
  switch(swtch) {
    case 1:
      if (swi_fp_trace_c != NULL) {
	fputs("\n[SWE_CALC]\n", swi_fp_trace_c);
	fprintf(swi_fp_trace_c, "  tjd = %.9f;", tjd);
	fprintf(swi_fp_trace_c, " ipl = %d;", ipl);
	fprintf(swi_fp_trace_c, " iflag = %d;\n", iflag);
	fprintf(swi_fp_trace_c, "  iflgret = swe_calc(tjd, ipl, iflag, xx, serr);");
	fprintf(swi_fp_trace_c, "	[ xx = %p ]\n", xx);
	fflush(swi_fp_trace_c);
      }
      break;
    case 2:
      if (swi_fp_trace_c != NULL) {
	fputs("  printf(\"swe_calc: %f\\t%d\\t%d\\t%f\\t%f\\t%f\\t%f\\t%f\\t%f\\t\", ", swi_fp_trace_c);
	fputs("\n\ttjd, ipl, iflgret, xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);\n", swi_fp_trace_c);
	fputs("  if (*serr != '\\0')", swi_fp_trace_c);
	fputs(" printf(serr);", swi_fp_trace_c);
	fputs(" printf(\"\\n\");\n", swi_fp_trace_c);
	fflush(swi_fp_trace_c);
      }
      if (swi_fp_trace_out != NULL) {
	fprintf(swi_fp_trace_out, "swe_calc: %f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t",
		      tjd, ipl, iflag, xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
	if (serr != NULL && *serr != '\0') {
	  fputs(serr, swi_fp_trace_out);
	}
	fputs("\n", swi_fp_trace_out);
	fflush(swi_fp_trace_out);
      }
      break;
    default:
      break;
  }
}

static void trace_swe_fixstar(int swtch, char *star, double tjd, int32 iflag, double *xx, char *serr)
{
  if (swi_trace_count >= TRACE_COUNT_MAX)
    return;
  switch(swtch) {
  case 1:
    if (swi_fp_trace_c != NULL) {
      fputs("\n[SWE_FIXSTAR]\n", swi_fp_trace_c);
      fprintf(swi_fp_trace_c, "  strcpy(star, \"%s\");", star);
      fprintf(swi_fp_trace_c, " tjd = %.9f;", tjd);
      fprintf(swi_fp_trace_c, " iflag = %d;\n", iflag);
      fprintf(swi_fp_trace_c, "  iflgret = swe_fixstar(star, tjd, iflag, xx, serr);");
      fprintf(swi_fp_trace_c, "   [ xx = %p ]\n", xx);
      fflush(swi_fp_trace_c);
    }
    break;
  case 2:
    if (swi_fp_trace_c != NULL) {
      fputs("  printf(\"swe_fixstar: %s\\t%f\\t%d\\t%f\\t%f\\t%f\\t%f\\t%f\\t%f\\t\", ", swi_fp_trace_c);
      fputs("\n\tstar, tjd, iflgret, xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);\n", swi_fp_trace_c);
      fputs("  if (*serr != '\\0')", swi_fp_trace_c);
      fputs(" printf(serr);", swi_fp_trace_c);
      fputs(" printf(\"\\n\");\n", swi_fp_trace_c);
      fflush(swi_fp_trace_c);
    }
    if (swi_fp_trace_out != NULL) {
      fprintf(swi_fp_trace_out, "swe_fixstar: %s\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t",
		    star, tjd, iflag, xx[0], xx[1], xx[2], xx[3], xx[4], xx[5]);
      if (serr != NULL && *serr != '\0') {
	fputs(serr, swi_fp_trace_out);
      }
      fputs("\n", swi_fp_trace_out);
      fflush(swi_fp_trace_out);
    }
    break;
  default:
    break;
  }
}

static void trace_swe_get_planet_name(int swtch, int ipl, char *s)
{
  if (swi_trace_count >= TRACE_COUNT_MAX)
    return;
  switch(swtch) {
    case 1:
      if (swi_fp_trace_c != NULL) {
	fputs("\n[SWE_GET_PLANET_NAME]\n", swi_fp_trace_c);
	fprintf(swi_fp_trace_c, "  ipl = %d;\n", ipl);
	fprintf(swi_fp_trace_c, "  swe_get_planet_name(ipl, s);");
	fprintf(swi_fp_trace_c, "   [ s = %p ]\n", s);
	fflush(swi_fp_trace_c);
      }
      break;
    case 2:
      if (swi_fp_trace_c != NULL) {
	fputs("  printf(\"swe_get_planet_name: %d\\t%s\\t\\n\", ", swi_fp_trace_c);
	fputs("ipl, s);\n", swi_fp_trace_c);
	fflush(swi_fp_trace_c);
      }
      if (swi_fp_trace_out != NULL) {
	fprintf(swi_fp_trace_out, "swe_get_planet_name: %d\t%s\t\n", ipl, s);
	fflush(swi_fp_trace_out);
      }
      break;
    default:
      break;
  }
}

#endif

void CALL_CONV swe_set_topo(double geolon, double geolat, double geoalt)
{
  swi_init_swed_if_start();
  if (swed.geopos_is_set == TRUE
    && swed.topd.geolon == geolon
    && swed.topd.geolat == geolat
    && swed.topd.geoalt == geoalt) {
    return;
  }
  swed.topd.geolon = geolon;
  swed.topd.geolat = geolat;
  swed.topd.geoalt = geoalt;
  swed.geopos_is_set = TRUE;

  swed.topd.teval = 0;

  swi_force_app_pos_etc();
}

void swi_force_app_pos_etc(void)
{
  int i;
  for (i = 0; i < SEI_NPLANETS; i++)
    swed.pldat[i].xflgs = -1;
  for (i = 0; i < SEI_NNODE_ETC; i++)
    swed.nddat[i].xflgs = -1;
  for (i = 0; i <= SE_NPLANETS; i++) {
    swed.savedat[i].tsave = 0;
    swed.savedat[i].iflgsave = -1;
  }
}

int swi_get_observer(double tjd, int32 iflag,
	AS_BOOL do_save, double *xobs, char *serr)
{
  int i;
  double sidt, delt, tjd_ut, eps, nut, nutlo[2];
  double f = EARTH_OBLATENESS;
  double re = EARTH_RADIUS;
  double cosfi, sinfi, cc, ss, cosl, sinl, h;
  if (!swed.geopos_is_set) {
    if (serr != NULL)
      strcpy(serr, "geographic position has not been set");
    return ERR;
  }

  delt = swe_deltat_ex(tjd, iflag, serr);
  tjd_ut = tjd - delt;
  if (swed.oec.teps == tjd && swed.nut.tnut == tjd) {
    eps = swed.oec.eps;
    nutlo[1] = swed.nut.nutlo[1];
    nutlo[0] = swed.nut.nutlo[0];
  } else {
    eps = swi_epsiln(tjd, iflag);
    if (!(iflag & SEFLG_NONUT))
      swi_nutation(tjd, iflag, nutlo);
  }
  if (iflag & SEFLG_NONUT) {
    nut = 0;
  } else {
    eps += nutlo[1];
    nut = nutlo[0];
  }

  sidt = swe_sidtime0(tjd_ut, eps * RADTODEG, nut * RADTODEG);
  sidt *= 15;

  cosfi = cos(swed.topd.geolat * DEGTORAD);
  sinfi = sin(swed.topd.geolat * DEGTORAD);
  cc= 1 / sqrt(cosfi * cosfi + (1-f) * (1-f) * sinfi * sinfi);
  ss= (1-f) * (1-f) * cc;

  cosl = cos((swed.topd.geolon + sidt) * DEGTORAD);
  sinl = sin((swed.topd.geolon + sidt) * DEGTORAD);
  h = swed.topd.geoalt;
  xobs[0] = (re * cc + h) * cosfi * cosl;
  xobs[1] = (re * cc + h) * cosfi * sinl;
  xobs[2] = (re * ss + h) * sinfi;

  swi_cartpol(xobs, xobs);

  xobs[3] = EARTH_ROT_SPEED;
  xobs[4] = xobs[5] = 0;
  swi_polcart_sp(xobs, xobs);

  for (i = 0; i <= 5; i++)
    xobs[i] /= AUNIT;

  if (!(iflag & SEFLG_NONUT)) {
    swi_coortrf2(xobs, xobs, -swed.nut.snut, swed.nut.cnut);

      swi_coortrf2(xobs+3, xobs+3, -swed.nut.snut, swed.nut.cnut);
    swi_nutate(xobs, iflag | SEFLG_SPEED, TRUE);
  }

  swi_precess(xobs, tjd, iflag, J_TO_J2000);

    swi_precess_speed(xobs, tjd, iflag, J_TO_J2000);

  if (do_save) {
    for (i = 0; i <= 5; i++)
      swed.topd.xobs[i] = xobs[i];
    swed.topd.teval = tjd;
    swed.topd.tjd_ut = tjd_ut;
  }
  return OK;
}

int32 CALL_CONV swe_time_equ(double tjd_ut, double *E, char *serr)
{
  int32 retval;
  double t, dt, x[6];
  double sidt = swe_sidtime(tjd_ut);
  int32 iflag = SEFLG_EQUATORIAL;
  iflag = plaus_iflag(iflag, -1, tjd_ut, serr);
  if (swi_init_swed_if_start() == 1 && !(iflag & SEFLG_MOSEPH) && serr != NULL) {
    strcpy(serr, "Please call swe_set_ephe_path() or swe_set_jplfile() before calling swe_time_equ(), swe_lmt_to_lat() or swe_lat_to_lmt()");
  }
  if (swed.jpl_file_is_open)
    iflag |= SEFLG_JPLEPH;
  t = tjd_ut + 0.5;
  dt = t - floor(t);
  sidt -= dt * 24;
  sidt *= 15;
  if ((retval = swe_calc_ut(tjd_ut, SE_SUN, iflag, x, serr)) == ERR) {
    *E = 0;
    return ERR;
  }
  dt = swe_degnorm(sidt - x[0] - 180);
  if (dt > 180)
    dt -= 360;
  dt *= 4;
  *E = dt / 1440.0;
  return OK;
}

int32 CALL_CONV swe_lmt_to_lat(double tjd_lmt, double geolon, double *tjd_lat, char *serr)
{
  int32 retval;
  double E, tjd_lmt0;
  tjd_lmt0 = tjd_lmt - geolon / 360.0;
  retval = swe_time_equ(tjd_lmt0, &E, serr);
  *tjd_lat = tjd_lmt + E;
  return retval;
}

int32 CALL_CONV swe_lat_to_lmt(double tjd_lat, double geolon, double *tjd_lmt, char *serr)
{
  int32 retval;
  double E, tjd_lmt0;
  tjd_lmt0 = tjd_lat - geolon / 360.0;
  retval = swe_time_equ(tjd_lmt0, &E, serr);

  retval = swe_time_equ(tjd_lmt0 - E, &E, serr);
  retval = swe_time_equ(tjd_lmt0 - E, &E, serr);
  *tjd_lmt = tjd_lat - E;
  return retval;
}

static int open_jpl_file(double *ss, char *fname, char *fpath, char *serr)
{
  int retc;
  char serr2[AS_MAXCH];
  retc = swi_open_jpl_file(ss, fname, fpath, serr);

  if (retc != OK && strstr(fname, SE_FNAME_DFT) != NULL && serr != NULL) {
    retc = swi_open_jpl_file(ss, SE_FNAME_DFT2, fpath, serr2);
    if (retc == OK) {
      strcpy(swed.jplfnam, SE_FNAME_DFT2);
      if (serr != NULL) {
        strcpy(serr2, "Error with JPL ephemeris file ");
	if (strlen(serr2) + strlen(SE_FNAME_DFT) < AS_MAXCH)
	  strcat(serr2, SE_FNAME_DFT);
	if (strlen(serr2) + strlen(serr) + 2 < AS_MAXCH)
	  sprintf(serr2 + strlen(serr2), ": %s", serr);
	if (strlen(serr2) + 17 < AS_MAXCH)
	  strcat(serr2, ". Defaulting to ");
	if (strlen(serr2) + strlen(SE_FNAME_DFT2) < AS_MAXCH)
	  strcat(serr2, SE_FNAME_DFT2);
        strcpy(serr, serr2);
      }
    }
  }
  if (retc == OK) {
    swed.jpldenum = swi_get_jpl_denum();
    swed.jpl_file_is_open = TRUE;
    swi_set_tid_acc(0, 0, swed.jpldenum, serr);
  }
  return retc;
}

#if 1
static int32 swi_fixstar_load_record(char *star, char *srecord, char *sname, char *sbayer, double *dparams, char *serr)
{
  char s[AS_MAXCH + 20], *sp, *sp2;
  char sstar[SWI_STAR_LENGTH + 1];
  char fstar[SWI_STAR_LENGTH + 1];
  int i, star_nr = 0;
  int line = 0;
  int fline = 0;
  int32 retc = OK;
  AS_BOOL  is_bayer = FALSE;
  size_t cmplen;
  size_t slen;
  struct fixed_star stardata;

  retc = fixstar_format_search_name(star, sstar, serr);
  if (retc == ERR)
    return ERR;

  if (*sstar == ',') {
    is_bayer = TRUE;

  } else if (isdigit((int) *sstar)) {
    star_nr = atoi(sstar);

  } else {
    if ((sp = strchr(sstar, ',')) != NULL)
      *sp = '\0';
  }
  cmplen = strlen(sstar);

  if (swed.fixfp == NULL) {
    if ((swed.fixfp = swi_fopen(SEI_FILE_FIXSTAR, SE_STARFILE, swed.ephepath, serr)) == NULL) {
      swed.is_old_starfile = TRUE;
      if ((swed.fixfp = swi_fopen(SEI_FILE_FIXSTAR, SE_STARFILE_OLD, swed.ephepath, NULL)) == NULL) {
	swed.is_old_starfile = FALSE;

	return ERR;
      }
    }
  }
  rewind(swed.fixfp);
  while (fgets(s, AS_MAXCH, swed.fixfp) != NULL) {
    fline++;

    if (*s == '#') continue;
    line++;

    if (star_nr == line) {
      goto found;
    } else if (star_nr > 0) {
      continue;
    }

    if ((sp = strchr(s, ',')) == NULL) {
      if (serr != NULL) {
	sprintf(serr, "star file %s damaged at line %d", SE_STARFILE, fline);
      }
      return ERR;
    }

    if (is_bayer) {
      if (strncmp(sp, sstar, cmplen) == 0) {
        goto found;
      } else {
        continue;
      }
    }

    *sp = '\0';

    slen = strlen(s);
    if (slen > SE_MAX_STNAME) slen = SE_MAX_STNAME;
    memcpy(fstar, s, slen);
    fstar[slen] = '\0';
    *sp = ',';

    while ((sp = strchr(fstar, ' ')) != NULL)
      swi_strcpy(sp, sp+1);
    i = (int) strlen(fstar);

    if (i < (int) cmplen)
      continue;

    for (sp2 = fstar; *sp2 != '\0'; sp2++) {
      *sp2 = tolower((int) *sp2);
    }
    if (strncmp(fstar, sstar, cmplen) == 0)
      goto found;
  }
  if (serr != NULL) {
    sprintf(serr, "star  not found");
    if (strlen(serr) + strlen(star) < AS_MAXCH) {
      sprintf(serr, "star %s not found", star);
    }
    return ERR;
  }
  found:
  strcpy(srecord, s);
  retc = fixstar_cut_string(srecord, star, &stardata, serr);
  if (retc == ERR) return ERR;
  if (dparams != NULL) {
    dparams[0] = stardata.epoch;

    dparams[1] = stardata.ra;

    dparams[2] = stardata.de;

    dparams[3] = stardata.ramot;

    dparams[4] = stardata.demot;

    dparams[5] = stardata.radvel;

    dparams[6] = stardata.parall;

    dparams[7] = stardata.mag;
  }
  return OK;
}

static int32 swi_fixstar_calc_from_record(char *srecord, double tjd, int32 iflag, char *star, double *xx, char *serr)
{
  int i;
  int32 retc = OK;
  double epoch, radv, parall;
  double ra_pm, de_pm, ra, de, t;
  struct fixed_star stardata;
  double daya, rdist;
  double x[6], xxsv[6], xobs[6], xobs_dt[6], *xpo = NULL, *xpo_dt = NULL;
  static TLS double xearth[6], xearth_dt[6], xsun[6], xsun_dt[6];
  double dt = PLAN_SPEED_INTV * 0.1;
  int32 epheflag, iflgsave;

  struct epsilon *oe = &swed.oec2000;
  iflgsave = iflag;
  iflag |= SEFLG_SPEED;
  if (serr != NULL)
    *serr = '\0';
  iflag = plaus_iflag(iflag, -1, tjd, serr);
  epheflag = iflag & SEFLG_EPHMASK;
  if (swi_init_swed_if_start() == 1 && !(epheflag & SEFLG_MOSEPH) && serr != NULL) {
    strcpy(serr, "Please call swe_set_ephe_path() or swe_set_jplfile() before calling swe_fixstar() or swe_fixstar_ut()");
  }
  if (swed.last_epheflag != epheflag) {
    free_planets();

    if (swed.jpl_file_is_open) {
      swi_close_jpl_file();
      swed.jpl_file_is_open = FALSE;
    }
    for (i = 0; i < SEI_NEPHFILES; i ++) {
      if (swed.fidat[i].fptr != NULL)
	fclose(swed.fidat[i].fptr);
      memset((void *) &swed.fidat[i], 0, sizeof(struct file_data));
    }
    swed.last_epheflag = epheflag;
  }

  if (iflag & SEFLG_SIDEREAL && !swed.ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);

  swi_check_ecliptic(tjd, iflag);

  swi_check_nutation(tjd, iflag);
  retc = fixstar_cut_string(srecord, star, &stardata, serr);
  if (retc == ERR) return ERR;
  epoch = stardata.epoch;
  ra_pm = stardata.ramot; de_pm = stardata.demot;
  radv = stardata.radvel; parall = stardata.parall;
  ra = stardata.ra; de = stardata.de;
  if (epoch == 1950) {
    t= (tjd - B1950);
  } else {
    t= (tjd - J2000);
  }
  x[0] = ra;
  x[1] = de;
  x[2] = 1;
  if (parall == 0) {
    rdist = 1000000000;
  } else {
    rdist = 1.0 / (parall * RADTODEG * 3600) * PARSEC_TO_AUNIT;

  }

  x[2] = rdist;
  x[3] = ra_pm / 36525.0;
  x[4] = de_pm / 36525.0;
  x[5] = radv / 36525.0;

  swi_polcart_sp(x, x);

  if (epoch == 1950) {
    swi_FK4_FK5(x, B1950);
    swi_precess(x, B1950, 0, J_TO_J2000);
    swi_precess(x+3, B1950, 0, J_TO_J2000);
  }

  if (epoch != 0) {
    swi_icrs2fk5(x, iflag, TRUE);

    if (swi_get_denum(SEI_SUN, iflag) >= 403) {
      swi_bias(x, J2000, SEFLG_SPEED, FALSE);
    }
  }

  if (!(iflag & SEFLG_BARYCTR) && (!(iflag & SEFLG_HELCTR) || !(iflag & SEFLG_MOSEPH))) {
    if ((retc =  main_planet_bary(tjd - dt, SEI_EARTH, epheflag, iflag, NO_SAVE, xearth_dt, xearth_dt, xsun_dt, NULL, serr)) != OK) {
      return ERR;
    }
    if ((retc =  main_planet_bary(tjd, SEI_EARTH, epheflag, iflag, DO_SAVE, xearth, xearth, xsun, NULL, serr)) != OK) {
      return ERR;
    }
  }

  if (iflag & SEFLG_TOPOCTR) {
    if (swi_get_observer(tjd - dt, iflag | SEFLG_NONUT, NO_SAVE, xobs_dt, serr) != OK)
      return ERR;
    if (swi_get_observer(tjd, iflag | SEFLG_NONUT, NO_SAVE, xobs, serr) != OK)
      return ERR;

    for (i = 0; i <= 5; i++) {
      xobs[i] = xobs[i] + xearth[i];
      xobs_dt[i] = xobs_dt[i] + xearth_dt[i];
    }
  } else if (!(iflag & SEFLG_BARYCTR) && (!(iflag & SEFLG_HELCTR) || !(iflag & SEFLG_MOSEPH))) {

    for (i = 0; i <= 5; i++) {
      xobs[i] = xearth[i];
      xobs_dt[i] = xearth_dt[i];
    }
  }

  if ((iflag & SEFLG_HELCTR) && (iflag & SEFLG_MOSEPH)) {
    xpo = NULL;
    xpo_dt = NULL;
  } else if (iflag & SEFLG_HELCTR) {
    xpo = xsun;
    xpo_dt = xsun_dt;
  } else if (iflag & SEFLG_BARYCTR) {
    xpo = NULL;
    xpo_dt = NULL;
  } else {
    xpo = xobs;
    xpo_dt = xobs_dt;
  }
  if (xpo == NULL) {
    for (i = 0; i <= 2; i++) {
      x[i] += t * x[i+3];
    }
  } else {
    for (i = 0; i <= 2; i++) {
      x[i] += t * x[i+3];
      x[i] -= xpo[i];
      x[i+3] -= xpo[i+3];
    }

  }

  if ((iflag & SEFLG_TRUEPOS) == 0 && (iflag & SEFLG_NOGDEFL) == 0) {
    swi_deflect_light(x, 0, iflag & SEFLG_SPEED);
  }

  if ((iflag & SEFLG_TRUEPOS) == 0 && (iflag & SEFLG_NOABERR) == 0)
    swi_aberr_light_ex(x, xpo, xpo_dt, dt, iflag & SEFLG_SPEED);

  if (!(iflag & SEFLG_ICRS) && (swi_get_denum(SEI_SUN, iflag) >= 403 || (iflag & SEFLG_BARYCTR))) {
    swi_bias(x, tjd, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = x[i];

  if ((iflag & SEFLG_J2000) == 0) {
    swi_precess(x, tjd, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED) {
      swi_precess_speed(x, tjd, iflag, J2000_TO_J);
    }
    oe = &swed.oec;
  } else {
    oe = &swed.oec2000;
  }

  if (!(iflag & SEFLG_NONUT))
    swi_nutate(x, iflag, FALSE);
if ((0)) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  printf("%.17f %.17f %f\n", x[0]/r, x[1]/r, x[2]/r);
}

  if ((iflag & SEFLG_EQUATORIAL) == 0) {
    swi_coortrf2(x, x, oe->seps, oe->ceps);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(x+3, x+3, oe->seps, oe->ceps);
    if (!(iflag & SEFLG_NONUT)) {
      swi_coortrf2(x, x, swed.nut.snut, swed.nut.cnut);
      if (iflag & SEFLG_SPEED)
	swi_coortrf2(x+3, x+3, swed.nut.snut, swed.nut.cnut);
    }
  }

  if (iflag & SEFLG_SIDEREAL) {

    if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
      if (swi_trop_ra2sid_lon(xxsv, x, xxsv, iflag) != OK)
        return ERR;
      if (iflag & SEFLG_EQUATORIAL) {
        for (i = 0; i <= 5; i++)
          x[i] = xxsv[i];
      }

    } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
      if (swi_trop_ra2sid_lon_sosy(xxsv, x, iflag) != OK)
	return ERR;
      if (iflag & SEFLG_EQUATORIAL) {
        for (i = 0; i <= 5; i++)
          x[i] = xxsv[i];
      }

    } else {
      swi_cartpol_sp(x, x);
      if (swi_get_ayanamsa_ex(tjd, iflag, &daya, serr) == ERR)
        return ERR;
      x[0] -= daya * DEGTORAD;
      swi_polcart_sp(x, x);
    }
  }

  if ((iflag & SEFLG_XYZ) == 0)
    swi_cartpol_sp(x, x);

  if ((iflag & SEFLG_RADIANS) == 0 && (iflag & SEFLG_XYZ) == 0) {
    for (i = 0; i < 2; i++) {
      x[i] *= RADTODEG;
      x[i+3] *= RADTODEG;
    }
  }
  for (i = 0; i <= 5; i++)
    xx[i] = x[i];
  if (!(iflgsave & SEFLG_SPEED)) {
    iflag = iflag & ~SEFLG_SPEED;
    for (i = 3; i <= 5; i++)
      xx[i] = 0;
  }

  if ((iflgsave & SEFLG_EPHMASK) == 0)
    iflag = iflag & ~SEFLG_DEFAULTEPH;
  return iflag;
}

int32 CALL_CONV swe_fixstar(char *star, double tjd, int32 iflag,
  double *xx, char *serr)
{
  int i;
  char sstar[SWI_STAR_LENGTH + 1];
  static TLS char slast_stardata[AS_MAXCH];
  static TLS char slast_starname[AS_MAXCH];
  char srecord[AS_MAXCH + 20], *sp;
  int retc;
  if (serr != NULL)
    *serr = '\0';
#ifdef TRACE
  swi_open_trace(serr);
  trace_swe_fixstar(1, star, tjd, iflag, xx, serr);
#endif
  retc = fixstar_format_search_name(star, sstar, serr);
  if (retc == ERR)
    goto return_err;
  if (*sstar == ',') {
    ;
  } else if (isdigit((int) *sstar)) {
    ;
  } else {
    if ((sp = strchr(sstar, ',')) != NULL)
      *sp = '\0';
  }

  if (*slast_stardata != '\0' && strcmp(slast_starname, sstar) == 0) {
    strcpy(srecord, slast_stardata);
    goto found;
  }
  if (get_builtin_star(star, sstar, srecord)) {
    goto found;
  }

  if ((retc = swi_fixstar_load_record(star, srecord, NULL, NULL, NULL, serr)) != OK)
    goto return_err;
  found:
  strcpy(slast_stardata, srecord);
  strcpy(slast_starname, sstar);
  if ((retc = swi_fixstar_calc_from_record(srecord, tjd, iflag, star, xx, serr)) == ERR)
    goto return_err;
#ifdef TRACE
  trace_swe_fixstar(2, star, tjd, iflag, xx, serr);
#endif
  return iflag;
  return_err:
  for (i = 0; i <= 5; i++)
    xx[i] = 0;
#ifdef TRACE
  trace_swe_fixstar(2, star, tjd, iflag, xx, serr);
#endif
  return retc;
}

int32 CALL_CONV swe_fixstar_ut(char *star, double tjd_ut, int32 iflag,
  double *xx, char *serr)
{
  double deltat;
  int32 retflag;
  int32 epheflag = 0;
  iflag = plaus_iflag(iflag, -1, tjd_ut, serr);
  epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag == 0) {
    epheflag = SEFLG_SWIEPH;
    iflag |= SEFLG_SWIEPH;
  }
  deltat = swe_deltat_ex(tjd_ut, iflag, serr);

  retflag = swe_fixstar(star, tjd_ut + deltat, iflag, xx, serr);
  if (retflag != ERR && (retflag & SEFLG_EPHMASK) != epheflag) {
    deltat = swe_deltat_ex(tjd_ut, retflag, NULL);
    retflag = swe_fixstar(star, tjd_ut + deltat, iflag, xx, NULL);
  }
  return retflag;
}

int32 CALL_CONV swe_fixstar_mag(char *star, double *mag, char *serr)
{
  char sstar[SWI_STAR_LENGTH + 1];
  static TLS char slast_stardata[AS_MAXCH];
  static TLS char slast_starname[AS_MAXCH];
  char srecord[AS_MAXCH + 20], *sp;
  struct fixed_star stardata;
  int retc;
  double dparams[20];
  if (serr != NULL)
    *serr = '\0';
  retc = fixstar_format_search_name(star, sstar, serr);
  if (retc == ERR)
    goto return_err;
  if (*sstar == ',') {
    ;
  } else if (isdigit((int) *sstar)) {
    ;
  } else {
    if ((sp = strchr(sstar, ',')) != NULL)
      *sp = '\0';
  }

  if (*slast_stardata != '\0' && strcmp(slast_starname, sstar) == 0) {
    strcpy(srecord, slast_stardata);
    retc = fixstar_cut_string(srecord, star, &stardata, serr);
    if (retc == ERR) goto return_err;

    dparams[7] = stardata.mag;
    goto found;
  }

  if ((retc = swi_fixstar_load_record(star, srecord, NULL, NULL, dparams, serr)) != OK)
    goto return_err;
  found:
  strcpy(slast_stardata, srecord);
  strcpy(slast_starname, sstar);
  *mag = dparams[7];
  return OK;
  return_err:
  *mag = 0;
  return retc;
}

#endif

int32 CALL_CONV swe_calc_pctr(double tjd, int32 ipl, int32 iplctr, int32 iflag, double *xxret, char *serr)
{
  double t = 0, dt, daya[2], dtsave_for_defl = 0;
  double xx[6], xxctr[6], xxctr2[6], xx0[6], xxsv[24], xxsp[6], dx[6], xreturn[24];
  double *xs;
  int i, j, niter;
  int32 iflag2, epheflag, retc;
  struct epsilon *oe;
  if (ipl == iplctr) {
    if (serr != NULL)
	  sprintf(serr, "ipl and iplctr (= %d) must not be identical\n", ipl);
	return ERR;
  }
  iflag = plaus_iflag(iflag, ipl, tjd, serr);
  epheflag = iflag & SEFLG_EPHMASK;

  swe_calc(tjd + swe_deltat_ex(tjd, epheflag, serr), SE_ECL_NUT, iflag, xx, serr);
  iflag &= ~(SEFLG_HELCTR|SEFLG_BARYCTR);
  iflag2 = epheflag;
  iflag2 |= (SEFLG_BARYCTR|SEFLG_J2000|SEFLG_ICRS|SEFLG_TRUEPOS|SEFLG_EQUATORIAL|SEFLG_XYZ|SEFLG_SPEED);
  iflag2 |= (SEFLG_NOABERR|SEFLG_NOGDEFL);
  retc = swe_calc(tjd, iplctr, iflag2, xxctr, serr);
  if (retc == ERR)
    return ERR;
  retc = swe_calc(tjd, ipl, iflag2, xx, serr);
  if (retc == ERR)
    return ERR;
  for (i = 0; i <= 5; i++) {
    xx0[i] = xx[i];

  }

  if (!(iflag & SEFLG_TRUEPOS)) {

    niter = 1;
    if (iflag & SEFLG_SPEED) {

      for (i = 0; i <= 2; i++)
	    xxsv[i] = xxsp[i] = xx[i] - xx[i+3];
      for (j = 0; j <= niter; j++) {
        for (i = 0; i <= 2; i++) {
          dx[i] = xxsp[i];
          dx[i] -= (xxctr[i] - xxctr[i+3]);
        }

        dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
        for (i = 0; i <= 2; i++)
          xxsp[i] = xxsv[i] - dt * xx0[i+3];
      }

      for (i = 0; i <= 2; i++)
        xxsp[i] = xxsv[i] - xxsp[i];
    }

    for (j = 0; j <= niter; j++) {
      for (i = 0; i <= 2; i++) {
        dx[i] = xx[i];
    	dx[i] -= xxctr[i];
      }
      dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;

      t = tjd - dt;
      dtsave_for_defl = dt;
      for (i = 0; i <= 2; i++)
        xx[i] = xx0[i] - dt * xx0[i+3];
    }

    if (iflag & SEFLG_SPEED) {
      for (i = 0; i <= 2; i++)
        xxsp[i] = xx0[i] - xx[i] - xxsp[i];
    }
    retc = swe_calc(t, iplctr, iflag2, xxctr2, serr);
    retc = swe_calc(t, ipl, iflag2, xx, serr);
  }

  if (!(iflag & SEFLG_HELCTR) && !(iflag & SEFLG_BARYCTR)) {

    for (i = 0; i <= 5; i++)
      xx[i] -= xxctr[i];
    if ((iflag & SEFLG_TRUEPOS) == 0 ) {

      if (iflag & SEFLG_SPEED)
        for (i = 3; i <= 5; i++)
          xx[i] -= xxsp[i-3];
    }
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOGDEFL))

    swi_deflect_light(xx, dtsave_for_defl, iflag);

  if (!(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOABERR)) {

    swi_aberr_light(xx, xxctr, iflag);

    if (iflag & SEFLG_SPEED) {
      for (i = 3; i <= 5; i++)
        xx[i] += xxctr[i] - xxctr2[i];
    }
  }
  if (!(iflag & SEFLG_SPEED))
    for (i = 3; i <= 5; i++)
      xx[i] = 0;

  if (!(iflag & SEFLG_ICRS) && swi_get_denum(ipl, epheflag) >= 403) {
    swi_bias(xx, t, iflag, FALSE);
  }

  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];

  if (!(iflag & SEFLG_J2000)) {
    swi_precess(xx, tjd, iflag, J2000_TO_J);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xx, tjd, iflag, J2000_TO_J);
    oe = &swed.oec;
  } else {
    oe = &swed.oec2000;
  }

  if (!(iflag & SEFLG_NONUT))
    swi_nutate(xx, iflag, FALSE);

  for (i = 0; i <= 5; i++)
    xreturn[18+i] = xx[i];

  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  if (iflag & SEFLG_SPEED)
    swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
  if (!(iflag & SEFLG_NONUT)) {
    swi_coortrf2(xx, xx, swed.nut.snut, swed.nut.cnut);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(xx+3, xx+3, swed.nut.snut, swed.nut.cnut);
  }

  for (i = 0; i <= 5; i++)
    xreturn[6+i] = xx[i];

  if (iflag & SEFLG_SIDEREAL) {

    if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
      if (swi_trop_ra2sid_lon(xxsv, xreturn+6, xreturn+18, iflag) != OK)
	return ERR;

    } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
      if (swi_trop_ra2sid_lon_sosy(xxsv, xreturn+6, iflag) != OK)
        return ERR;
    } else {

      swi_cartpol_sp(xreturn+6, xreturn);

      for (i = 0; i < 24; i++)
        xxsv[i] = xreturn[i];
      if (swi_get_ayanamsa_with_speed(tjd, iflag, daya, serr) == ERR)
        return ERR;

      for (i = 0; i < 24; i++)
        xreturn[i] = xxsv[i];
      xreturn[0] -= daya[0] * DEGTORAD;
      xreturn[3] -= daya[1] * DEGTORAD;
      swi_polcart_sp(xreturn, xreturn+6);
    }
  }

  swi_cartpol_sp(xreturn+18, xreturn+12);
  swi_cartpol_sp(xreturn+6, xreturn);

  for (i = 0; i < 2; i++) {
    xreturn[i] *= RADTODEG;
    xreturn[i+3] *= RADTODEG;
    xreturn[i+12] *= RADTODEG;
    xreturn[i+15] *= RADTODEG;
  }

  if (iflag & SEFLG_EQUATORIAL) {
    xs = xreturn+12;
  } else {
    xs = xreturn;
  }
  if (iflag & SEFLG_XYZ)
    xs = xs+6;
  for (i = 0; i < 6; i++)
    xxret[i] = xs[i];
  if (!(iflag & SEFLG_SPEED)) {
    for (i = 3; i < 6; i++)
      xxret[i] = 0;
  }
  if (iflag & SEFLG_RADIANS) {
    for (i = 0; i < 2; i++)
      xxret[i] *= DEGTORAD;
    if (iflag & SEFLG_SPEED) {
      for (i = 3; i < 5; i++)
	xxret[i] *= DEGTORAD;
    }
  }
  if (retc == ERR)
    return ERR;
  return(iflag);
}

const char *CALL_CONV swe_get_current_file_data(int ifno, double *tfstart, double *tfend, int *denum)
{
  if (ifno < 0 || ifno > 4) return NULL;
  struct file_data *pfp = &swed.fidat[ifno];
  if (strlen(pfp->fnam) == 0) return NULL;
  *tfstart = pfp->tfstart;
  *tfend = pfp->tfend;
  *denum = pfp->sweph_denum;
  return pfp->fnam;
}

#define CROSS_PRECISION (1 / 3600000.0)

double CALL_CONV swe_solcross(double x2cross, double jd_et, int flag, char *serr)
{
  double x[6], xlp, dist;
  double jd;
  int ipl = SE_SUN;

  flag |= SEFLG_SPEED;
  if (swe_calc(jd_et, ipl, flag, x, serr) < 0)
    return jd_et - 1;
  xlp = 360.0 / 365.24;
  dist = swe_degnorm(x2cross - x[0]);
  jd = jd_et + dist / xlp;
  for(;;) {
    if (swe_calc(jd, ipl, flag, x, serr) < 0)
      return jd_et - 1;
    dist = swe_difdeg2n(x2cross, x[0]);
    jd += dist / x[3];
    if (fabs(dist) < CROSS_PRECISION) break;
  }
  return jd;
}

double CALL_CONV swe_solcross_ut(double x2cross, double jd_ut, int flag, char *serr)
{
  double x[6], xlp, dist;
  double jd;
  int ipl = SE_SUN;

  flag |= SEFLG_SPEED;
  if (swe_calc_ut(jd_ut, ipl, flag, x, serr) < 0)
    return jd_ut - 1;
  xlp = 360.0 / 365.24;
  dist = swe_degnorm(x2cross - x[0]);
  jd = jd_ut + dist / xlp;
  for(;;) {
    if (swe_calc_ut(jd, ipl, flag, x, serr) < 0)
      return jd_ut - 1;
    dist = swe_difdeg2n(x2cross, x[0]);
    jd += dist / x[3];
    if (fabs(dist) < CROSS_PRECISION) break;
  }
  return jd;
}

double CALL_CONV swe_mooncross(double x2cross, double jd_et, int flag, char *serr)
{
  double x[6], xlp, dist;
  double jd;
  int ipl = SE_MOON;

  flag |= SEFLG_SPEED;
  if (swe_calc(jd_et, ipl, flag, x, serr) < 0)
    return jd_et - 1;
  xlp = 360.0 / 27.32;
  dist = swe_degnorm(x2cross - x[0]);
  jd = jd_et + dist / xlp;
  for(;;) {
    if (swe_calc(jd, ipl, flag, x, serr) < 0)
      return jd_et - 1;
    dist = swe_difdeg2n(x2cross, x[0]);
    jd += dist / x[3];
    if (fabs(dist) < CROSS_PRECISION) break;
  }
  return jd;
}

double CALL_CONV swe_mooncross_ut(double x2cross, double jd_ut, int flag, char *serr)
{
  double x[6], xlp, dist;
  double jd;
  int ipl = SE_MOON;

  flag |= SEFLG_SPEED;
  if (swe_calc_ut(jd_ut, ipl, flag, x, serr) < 0)
    return jd_ut - 1;
  xlp = 360.0 / 27.32;
  dist = swe_degnorm(x2cross - x[0]);
  jd = jd_ut + dist / xlp;
  for(;;) {
    if (swe_calc_ut(jd, ipl, flag, x, serr) < 0)
      return jd_ut - 1;
    dist = swe_difdeg2n(x2cross, x[0]);
    jd += dist / x[3];
    if (fabs(dist) < CROSS_PRECISION) break;
  }
  return jd;
}

double CALL_CONV swe_mooncross_node(double jd_et, int flag, double *xlon, double *xla, char *serr)
{
  double x[6], xlat, dist;
  double jd;
  int ipl = SE_MOON;
  flag |= SEFLG_SPEED;
  if (swe_calc(jd_et, ipl, flag, x, serr) < 0)
    return jd_et - 1;
  xlat = x[1];
  jd = jd_et + 1;
  for(;;) {
    if (swe_calc(jd, ipl, flag, x, serr) < 0)
      return jd_et - 1;
    if ((x[1] >= 0 && xlat < 0) || (x[1] < 0 && xlat > 0))
      break;
    jd += 1;
  }
  dist = x[1];
  for(;;) {
    jd -= dist / x[4];
    if (swe_calc(jd, ipl, flag, x, serr) < 0)
      return jd_et - 1;
    dist = x[1];
    if (fabs(dist) < CROSS_PRECISION) {
      *xlon = x[0];
      *xla = x[1];
      break;
    }
  }
  return jd;
}

double CALL_CONV swe_mooncross_node_ut(double jd_ut, int flag, double *xlon, double *xla, char *serr)
{
  double x[6], xlat, dist;
  double jd;
  int ipl = SE_MOON;
  flag |= SEFLG_SPEED;
  if (swe_calc_ut(jd_ut, ipl, flag, x, serr) < 0)
    return jd_ut - 1;
  xlat = x[1];
  jd = jd_ut + 1;
  for(;;) {
    if (swe_calc_ut(jd, ipl, flag, x, serr) < 0)
      return jd_ut - 1;
    if ((x[1] >= 0 && xlat < 0) || (x[1] < 0 && xlat > 0))
      break;
    jd += 1;
  }
  dist = x[1];
  for(;;) {
    jd -= dist / x[4];
    if (swe_calc_ut(jd, ipl, flag, x, serr) < 0)
      return jd_ut - 1;
    dist = x[1];
    if (fabs(dist) < CROSS_PRECISION) {
      *xlon = x[0];
      *xla = x[1];
      break;
    }
  }
  return jd;
}

int32 CALL_CONV swe_helio_cross(int ipl, double x2cross, double jd_et, int iflag, int dir, double *jd_cross, char *serr)
{
  double x[6], xlp, dist;
  double jd;
  int flag = iflag | SEFLG_SPEED | SEFLG_HELCTR;
  if (ipl == SE_SUN
    || ipl == SE_MOON
    || (ipl >= SE_MEAN_NODE && ipl <= SE_OSCU_APOG)
    || (ipl >= SE_INTP_APOG && ipl < SE_NPLANETS)
  ) {
    char snam[AS_MAXCH];
    swe_get_planet_name(ipl, snam);
    if (serr != NULL) sprintf(serr, "swe_helio_cross: not possible for object %d = %s", ipl, snam);
    return ERR;
  }
  if (swe_calc(jd_et, ipl, flag, x, serr) < 0)
    return ERR;
  xlp = x[3];
  if (ipl == SE_CHIRON)
    xlp = 0.01971;
  dist = swe_degnorm(x2cross - x[0]);
  if (dir >= 0) {
    jd = jd_et + dist / xlp;
  } else {
    dist = 360.0 - dist;
    jd = jd_et - dist / xlp;
  }
  for(;;) {
    if (swe_calc(jd, ipl, flag, x, serr) < 0)
      return ERR;
    dist = swe_difdeg2n(x2cross, x[0]);
    jd += dist / x[3];
    if (fabs(dist) < CROSS_PRECISION) break;
  }
  *jd_cross = jd;
  return OK;
}

int32 CALL_CONV swe_helio_cross_ut(int ipl, double x2cross, double jd_ut, int iflag, int dir, double *jd_cross, char *serr)
{
  double x[6], xlp, dist;
  double jd;
  int flag = iflag | SEFLG_SPEED | SEFLG_HELCTR;
  if (ipl == SE_SUN
    || ipl == SE_MOON
    || (ipl >= SE_MEAN_NODE && ipl <= SE_OSCU_APOG)
    || (ipl >= SE_INTP_APOG && ipl < SE_NPLANETS)
  ) {
    char snam[AS_MAXCH];
    swe_get_planet_name(ipl, snam);
    if (serr != NULL) sprintf(serr, "swe_helio_cross: not possible for object %d = %s", ipl, snam);
    return ERR;
  }
  if (swe_calc_ut(jd_ut, ipl, flag, x, serr) < 0)
    return ERR;
  xlp = x[3];
  if (ipl == SE_CHIRON)
    xlp = 0.01971;
  dist = swe_degnorm(x2cross - x[0]);
  if (dir >= 0) {
    jd = jd_ut + dist / xlp;
  } else {
    dist = 360.0 - dist;
    jd = jd_ut - dist / xlp;
  }
  for(;;) {
    if (swe_calc_ut(jd, ipl, flag, x, serr) < 0)
      return ERR;
    dist = swe_difdeg2n(x2cross, x[0]);
    jd += dist / x[3];
    if (fabs(dist) < CROSS_PRECISION) break;
  }
  *jd_cross = jd;
  return OK;
}
