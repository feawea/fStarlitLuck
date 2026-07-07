#if MSDOS
#else
  #define _FILE_OFFSET_BITS 64
#endif

#include <string.h>
#include "swephexp.h"
#include "sweph.h"
#include "swejpl.h"

#if MSDOS
  typedef __int64 off_t64;
  #define FSEEK _fseeki64
  #define FTELL _ftelli64
#else
  typedef off_t off_t64;
  #define FSEEK fseeko
  #define FTELL ftello

#endif

#define DEBUG_DO_SHOW	FALSE

struct jpl_save {
  char *jplfname;
  char *jplfpath;
  FILE *jplfptr;
  short do_reorder;
  double eh_cval[400];
  double eh_ss[3], eh_au, eh_emrat;
  int32 eh_denum, eh_ncon, eh_ipt[39];
  char ch_cnam[6*400];
  double pv[78];
  double pvsun[6];
  double buf[1500];
  double pc[18], vc[18], ac[18], jc[18];
  short do_km;
};

static TLS struct jpl_save *js;

static int state (double et, int32 *list, int do_bary,
		  double *pv, double *pvsun, double *nut, char *serr);
static int interp(double *buf, double t, double intv, int32 ncfin,
		  int32 ncmin, int32 nain, int32 ifl, double *pv);
static int32 fsizer(char *serr);
static void reorder(char *x, int size, int number);
static int read_const_jpl(double *ss, char *serr);

static int32 fsizer(char *serr)
{

  int32 ncon;
  double emrat;
  int32 numde;
  double au, ss[3];
  int i, kmx, khi, nd;
  int32 ksize, lpt[3];
  char ttl[6*14*3];
  size_t nrd;
  if ((js->jplfptr = swi_fopen(SEI_FILE_PLANET, js->jplfname, js->jplfpath, serr)) == NULL) {
    return NOT_AVAILABLE;
  }

  nrd = fread((void *) &ttl[0], 1, 252, js->jplfptr);
  if (nrd != 252) return NOT_AVAILABLE;

  nrd = fread((void *) js->ch_cnam, 1, 6*400, js->jplfptr);
  if (nrd != 6*400) return NOT_AVAILABLE;

  nrd = fread((void *) &ss[0], sizeof(double), 3, js->jplfptr);
  if (nrd != 3) return NOT_AVAILABLE;

  if (ss[2] < 1 || ss[2] > 200)
    js->do_reorder = TRUE;
  else
    js->do_reorder = 0;
  for (i = 0; i < 3; i++)
    js->eh_ss[i] = ss[i];
  if (js->do_reorder)
    reorder((char *) &js->eh_ss[0], sizeof(double), 3);

  if (js->eh_ss[0] < -5583942 || js->eh_ss[1] > 9025909 || js->eh_ss[2] < 1 || js->eh_ss[2] > 200) {
    if (serr != NULL) {
      strcpy(serr, "alleged ephemeris file has invalid format.");
      if (strlen(serr) + strlen(js->jplfname) + 3 < AS_MAXCH) {
	sprintf(serr, "alleged ephemeris file (%s) has invalid format.", js->jplfname);
      }
    }
    return(NOT_AVAILABLE);
  }

  nrd = fread((void *) &ncon, sizeof(int32), 1, js->jplfptr);
  if (nrd != 1) return NOT_AVAILABLE;
  if (js->do_reorder)
    reorder((char *) &ncon, sizeof(int32), 1);

  nrd = fread((void *) &au, sizeof(double), 1, js->jplfptr);
  if (nrd != 1) return NOT_AVAILABLE;
  if (js->do_reorder)
    reorder((char *) &au, sizeof(double), 1);

  nrd = fread((void *) &emrat, sizeof(double), 1, js->jplfptr);
  if (nrd != 1) return NOT_AVAILABLE;
  if (js->do_reorder)
    reorder((char *) &emrat, sizeof(double), 1);

  nrd = fread((void *) &js->eh_ipt[0], sizeof(int32), 36, js->jplfptr);
  if (nrd != 36) return NOT_AVAILABLE;
  if (js->do_reorder)
    reorder((char *) &js->eh_ipt[0], sizeof(int32), 36);

  nrd = fread((void *) &numde, sizeof(int32), 1, js->jplfptr);
  if (nrd != 1) return NOT_AVAILABLE;
  if (js->do_reorder)
    reorder((char *) &numde, sizeof(int32), 1);

  nrd = fread(&lpt[0], sizeof(int32), 3, js->jplfptr);
  if (nrd != 3) return NOT_AVAILABLE;
  if (js->do_reorder)
    reorder((char *) &lpt[0], sizeof(int32), 3);

  for (i = 0; i < 3; ++i)
    js->eh_ipt[i + 36] = lpt[i];
  rewind(js->jplfptr);

  kmx = 0;
  khi = 0;
  for (i = 0; i < 13; i++) {
    if (js->eh_ipt[i * 3] > kmx) {
      kmx = js->eh_ipt[i * 3];
      khi = i + 1;
    }
  }
  if (khi == 12)
    nd = 2;
  else
    nd = 3;
  ksize = (js->eh_ipt[khi * 3 - 3] + nd * js->eh_ipt[khi * 3 - 2] * js->eh_ipt[khi * 3 - 1] - 1L) * 2L;

  if (ksize == 1546)
    ksize = 1652;
#if 0

  switch (numde) {
    case 403:
    case 405:
    case 410:
    case 413:
    case 414:
    case 418:
    case 421:
      ksize = 2036;
      break;
    case 404:
    case 406:
      ksize = 1456;
      break;
    case 200:
      ksize = 1652;
      break;
    case 102:
      ksize = 1652;
      break;
    default:
      if (serr != NULL)
	sprintf(serr,"unknown numde value %d;", numde);
      return ERR;
  }
#endif
  if (ksize < 1000 || ksize > 5000) {
    if (serr != NULL)
      sprintf(serr, "JPL ephemeris file does not provide valid ksize (%d)", ksize);
    return NOT_AVAILABLE;
  }
  return ksize;
}

int swi_pleph(double et, int ntarg, int ncent, double *rrd, char *serr)
{
  int i, retc;
  int32 list[12];
  double *pv = js->pv;
  double *pvsun = js->pvsun;
  for (i = 0; i < 6; ++i)
    rrd[i] = 0.0;
  if (ntarg == ncent)
      return 0;
  for (i = 0; i < 12; ++i)
    list[i] = 0;

  if (ntarg == J_NUT) {
    if (js->eh_ipt[34] > 0) {
      list[10] = 2;
      return(state(et, list, FALSE, pv, pvsun, rrd, serr));
    } else {
      if (serr != NULL)
	sprintf(serr,"No nutations on the JPL ephemeris file;");
      return (NOT_AVAILABLE);
    }
  }
  if (ntarg == J_LIB) {
    if (js->eh_ipt[37] > 0) {
      list[11] = 2;
      if ((retc = state(et, list, FALSE, pv, pvsun, rrd, serr)) != OK)
	return (retc);
      for (i = 0; i < 6; ++i)
	rrd[i] = pv[i + 60];
      return 0;
    } else {
      if (serr != NULL)
	sprintf(serr,"No librations on the ephemeris file;");
      return (NOT_AVAILABLE);
    }
  }

  if (ntarg < J_SUN)
    list[ntarg] = 2;
  if (ntarg == J_MOON)
    list[J_EARTH] = 2;
  if (ntarg == J_EARTH)
    list[J_MOON] = 2;
  if (ntarg == J_EMB)
    list[J_EARTH] = 2;
  if (ncent < J_SUN)
    list[ncent] = 2;
  if (ncent == J_MOON)
    list[J_EARTH] = 2;
  if (ncent == J_EARTH)
    list[J_MOON] = 2;
  if (ncent == J_EMB)
    list[J_EARTH] = 2;
  if ((retc = state(et, list, TRUE, pv, pvsun, rrd, serr)) != OK)
    return (retc);
  if (ntarg == J_SUN || ncent == J_SUN) {
    for (i = 0; i < 6; ++i)
      pv[i + 6*J_SUN] = pvsun[i];
  }
  if (ntarg == J_SBARY || ncent == J_SBARY) {
    for (i = 0; i < 6; ++i) {
      pv[i + 6*J_SBARY] = 0.;
    }
  }
  if (ntarg == J_EMB || ncent == J_EMB) {
    for (i = 0; i < 6; ++i)
      pv[i + 6*J_EMB] = pv[i + 6*J_EARTH];
  }
  if ((ntarg==J_EARTH && ncent==J_MOON) || (ntarg == J_MOON && ncent==J_EARTH)){
    for (i = 0; i < 6; ++i)
      pv[i + 6*J_EARTH] = 0.;

  } else {
    if (list[J_EARTH] == 2) {
      for (i = 0; i < 6; ++i)
	pv[i + 6*J_EARTH] -= pv[i + 6*J_MOON] / (js->eh_emrat + 1.);
    }
    if (list[J_MOON] == 2) {
      for (i = 0; i < 6; ++i) {
	pv[i + 6*J_MOON] += pv[i + 6*J_EARTH];
      }
    }
  }
  for (i = 0; i < 6; ++i)
    rrd[i] = pv[i + ntarg * 6] - pv[i + ncent * 6];
  return OK;
}

static int interp(double *buf, double t, double intv, int32 ncfin,
		  int32 ncmin, int32 nain, int32 ifl, double *pv)
{

  static TLS int np, nv;
  static TLS int nac;
  static TLS int njk;
  static TLS double twot = 0.;
  double *pc = js->pc;
  double *vc = js->vc;
  double *ac = js->ac;
  double *jc = js->jc;
  int ncf = (int) ncfin;
  int ncm = (int) ncmin;
  int na = (int) nain;

  double temp;
  int i, j, ni;
  double tc;
  double dt1, bma;
  double bma2, bma3;

  if (t >= 0)
    dt1 = floor(t);
  else
    dt1 = -floor(-t);
  temp = na * t;
  ni = (int) (temp - dt1);

  tc = (fmod(temp, 1.0) + dt1) * 2. - 1.;

  if (tc != pc[1]) {
    np = 2;
    nv = 3;
    nac = 4;
    njk = 5;
    pc[1] = tc;
    twot = tc + tc;
  }

  if (np < ncf) {
    for (i = np; i < ncf; ++i)
      pc[i] = twot * pc[i - 1] - pc[i - 2];
    np = ncf;
  }

  for (i = 0; i < ncm; ++i) {
    pv[i] = 0.;
    for (j = ncf-1; j >= 0; --j)
      pv[i] += pc[j] * buf[j + (i + ni * ncm) * ncf];
  }
  if (ifl <= 1)
    return 0;

  bma = (na + na) / intv;
  vc[2] = twot + twot;
  if (nv < ncf) {
    for (i = nv; i < ncf; ++i)
      vc[i] = twot * vc[i - 1] + pc[i - 1] + pc[i - 1] - vc[i - 2];
    nv = ncf;
  }

  for (i = 0; i < ncm; ++i) {
    pv[i + ncm] = 0.;
    for (j = ncf-1; j >= 1; --j)
      pv[i + ncm] += vc[j] * buf[j + (i + ni * ncm) * ncf];
    pv[i + ncm] *= bma;
  }
  if (ifl == 2)
    return 0;

  bma2 = bma * bma;
  ac[3] = pc[1] * 24.;
  if (nac < ncf) {
    nac = ncf;
    for (i = nac; i < ncf; ++i)
      ac[i] = twot * ac[i - 1] + vc[i - 1] * 4. - ac[i - 2];
  }

  for (i = 0; i < ncm; ++i) {
    pv[i + ncm * 2] = 0.;
    for (j = ncf-1; j >= 2; --j)
      pv[i + ncm * 2] += ac[j] * buf[j + (i + ni * ncm) * ncf];
    pv[i + ncm * 2] *= bma2;
  }
  if (ifl == 3)
      return 0;

  bma3 = bma * bma2;
  jc[4] = pc[1] * 192.;
  if (njk < ncf) {
    njk = ncf;
    for (i = njk; i < ncf; ++i)
      jc[i] = twot * jc[i - 1] + ac[i - 1] * 6. - jc[i - 2];
  }

  for (i = 0; i < ncm; ++i) {
    pv[i + ncm * 3] = 0.;
    for (j = ncf-1; j >= 3; --j)
      pv[i + ncm * 3] += jc[j] * buf[j + (i + ni * ncm) * ncf];
    pv[i + ncm * 3] *= bma3;
  }
  return 0;
}

static int state(double et, int32 *list, int do_bary,
	  double *pv, double *pvsun, double *nut, char *serr)
{
  int i, j, k;
  int32 nseg;
  off_t64 flen, nb;
  double *buf = js->buf;
  double aufac, s, t, intv, ts[4];
  int32 nrecl, ksize;
  int32 nr;
  double et_mn, et_fr;
  int32 *ipt = js->eh_ipt;
  char ch_ttl[252];
  static TLS int32 irecsz;
  static TLS int32 nrl, lpt[3], ncoeffs;
  size_t nrd;
  if (js->jplfptr == NULL) {
    ksize = fsizer(serr);
    nrecl = 4;
    if (ksize == NOT_AVAILABLE)
      return NOT_AVAILABLE;
    irecsz = nrecl * ksize;
    ncoeffs = ksize / 2;

    nrd = fread((void *) ch_ttl, 1, 252, js->jplfptr);
    if (nrd != 252) return NOT_AVAILABLE;

    nrd = fread((void *) js->ch_cnam, 1, 2400, js->jplfptr);
    if (nrd != 2400) return NOT_AVAILABLE;

    nrd = fread((void *) &js->eh_ss[0], sizeof(double), 3, js->jplfptr);
    if (nrd != 3) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &js->eh_ss[0], sizeof(double), 3);

    nrd = fread((void *) &js->eh_ncon, sizeof(int32), 1, js->jplfptr);
    if (nrd != 1) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &js->eh_ncon, sizeof(int32), 1);

    nrd = fread((void *) &js->eh_au, sizeof(double), 1, js->jplfptr);
    if (nrd != 1) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &js->eh_au, sizeof(double), 1);

    nrd = fread((void *) &js->eh_emrat, sizeof(double), 1, js->jplfptr);
    if (nrd != 1) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &js->eh_emrat, sizeof(double), 1);

    nrd = fread((void *) &ipt[0], sizeof(int32), 36, js->jplfptr);
    if (nrd != 36) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &ipt[0], sizeof(int32), 36);

    nrd = fread((void *) &js->eh_denum, sizeof(int32), 1, js->jplfptr);
    if (nrd != 1) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &js->eh_denum, sizeof(int32), 1);
    nrd = fread((void *) &lpt[0], sizeof(int32), 3, js->jplfptr);
    if (nrd != 3) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &lpt[0], sizeof(int32), 3);

    FSEEK(js->jplfptr, (off_t64) (1L * irecsz), 0);
    nrd = fread((void *) &js->eh_cval[0], sizeof(double), 400, js->jplfptr);
    if (nrd != 400) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &js->eh_cval[0], sizeof(double), 400);

    for (i = 0; i < 3; ++i)
      ipt[i + 36] = lpt[i];
    nrl = 0;

    FSEEK(js->jplfptr, (off_t64) 0L, SEEK_END);
    flen = FTELL(js->jplfptr);

    nseg = (int32) ((js->eh_ss[1] - js->eh_ss[0]) / js->eh_ss[2]);

    for(i = 0, nb = 0; i < 13; i++) {
      k = 3;
      if (i == 11)
	k = 2;
      nb += (ipt[i*3+1] * ipt[i*3+2]) * k * nseg;
    }

    nb += 2 * nseg;

    nb *= 8;

    nb += 2 * ksize * nrecl;
    if (flen != nb

      && flen - nb != ksize * nrecl
      ) {
      if (serr != NULL) {
	sprintf(serr, "JPL ephemeris file is mutilated; length = %d instead of %d.", (unsigned int) flen, (unsigned int) nb);
	if (strlen(serr) + strlen(js->jplfname) < AS_MAXCH - 1) {
	  sprintf(serr, "JPL ephemeris file %s is mutilated; length = %d instead of %d.", js->jplfname, (unsigned int) flen, (unsigned int) nb);
	}
      }
      return(NOT_AVAILABLE);
    }

    FSEEK(js->jplfptr, (off_t64) (2L * irecsz), 0);
    nrd = fread((void *) &ts[0], sizeof(double), 2, js->jplfptr);
    if (nrd != 2) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &ts[0], sizeof(double), 2);
    FSEEK(js->jplfptr, (off_t64) ((nseg + 2 - 1) * ((off_t64) irecsz)), 0);
    nrd = fread((void *) &ts[2], sizeof(double), 2, js->jplfptr);
    if (nrd != 2) return NOT_AVAILABLE;
    if (js->do_reorder)
      reorder((char *) &ts[2], sizeof(double), 2);
    if (ts[0] != js->eh_ss[0] || ts[3] != js->eh_ss[1]) {
      if (serr != NULL)
	sprintf(serr, "JPL ephemeris file is corrupt; start/end date check failed. %.1f != %.1f || %.1f != %.1f", ts[0],js->eh_ss[0],ts[3],js->eh_ss[1]);
      return NOT_AVAILABLE;
    }
  }
  if (list == NULL)
    return 0;
  s = et - .5;
  et_mn = floor(s);
  et_fr = s - et_mn;
  et_mn += .5;

  if (et < js->eh_ss[0] || et > js->eh_ss[1]) {
    if (serr != NULL)
      sprintf(serr,"jd %f outside JPL eph. range %.2f .. %.2f;", et, js->eh_ss[0], js->eh_ss[1]);
    return BEYOND_EPH_LIMITS;
  }

  nr = (int32) ((et_mn - js->eh_ss[0]) / js->eh_ss[2]) + 2;
  if (et_mn == js->eh_ss[1])
    --nr;
  t = (et_mn - ((nr - 2) * js->eh_ss[2] + js->eh_ss[0]) + et_fr) / js->eh_ss[2];

  if (nr != nrl) {
    nrl = nr;
    if (FSEEK(js->jplfptr, (off_t64) (nr * ((off_t64) irecsz)), 0) != 0) {
      if (serr != NULL)
	sprintf(serr, "Read error in JPL eph. at %f\n", et);
      return NOT_AVAILABLE;
    }
    for (k = 1; k <= ncoeffs; ++k) {
      if ( fread((void *) &buf[k - 1], sizeof(double), 1, js->jplfptr) != 1) {
	if (serr != NULL)
	  sprintf(serr, "Read error in JPL eph. at %f\n", et);
	return NOT_AVAILABLE;
      }
      if (js->do_reorder)
	reorder((char *) &buf[k-1], sizeof(double), 1);
    }
  }
  if (js->do_km) {
    intv = js->eh_ss[2] * 86400.;
    aufac = 1.;
  } else {
    intv = js->eh_ss[2];
    aufac = 1. / js->eh_au;
  }

  interp(&buf[(int) ipt[30] - 1], t, intv, ipt[31], 3L, ipt[32], 2L, pvsun);
  for (i = 0; i < 6; ++i) {
    pvsun[i] *= aufac;
  }

  for (i = 0; i < 10; ++i) {
    if (list[i] > 0) {
      interp(&buf[(int) ipt[i * 3] - 1], t, intv, ipt[i * 3 + 1], 3L,
	     ipt[i * 3 + 2], list[i], &pv[i * 6]);
      for (j = 0; j < 6; ++j) {
	if (i < 9 && ! do_bary) {
	  pv[j + i * 6] = pv[j + i * 6] * aufac - pvsun[j];
	} else {
	  pv[j + i * 6] *= aufac;
	}
      }
    }
  }

  if (list[10] > 0 && ipt[34] > 0) {
    interp(&buf[(int) ipt[33] - 1], t, intv, ipt[34], 2L, ipt[35],
	     list[10], nut);
  }

  if (list[11] > 0 && ipt[37] > 0) {
    interp(&buf[(int) ipt[36] - 1], t, intv, ipt[37], 3L, ipt[38], list[1],
	    &pv[60]);
  }
  return OK;
}

static int read_const_jpl(double *ss,  char *serr)
{
  int i, retc;
  retc = state(0.0, NULL, FALSE, NULL, NULL, NULL, serr);
  if (retc != OK)
    return (retc);
  for (i = 0; i < 3; i++)
    ss[i] = js->eh_ss[i];
#if DEBUG_DO_SHOW
  {
    static const char *bname[] = {
	"Mercury", "Venus", "EMB", "Mars", "Jupiter", "Saturn",
	"Uranus", "Neptune", "Pluto", "Moon", "SunBary", "Nut", "Libr"};
    int j, k;
    int32 nb, nc;
    printf(" JPL TEST-EPHEMERIS program.  Version October 1995.\n");
    for (i = 0; i < 13; i++) {
      j = i * 3;
      k = 3;
      if (i == 11) k = 2;
      nb = js->eh_ipt[j+1] * js->eh_ipt[j+2] * k;
      nc = (int32) (nb * 36525L / js->eh_ss[2] * 8L);
      printf("%s\t%d\tipt[%d]\t%3ld %2ld %2ld,\t",
	bname[i], i, j, js->eh_ipt[j], js->eh_ipt[j+1], js->eh_ipt[j+2]);
      printf("%3ld double, bytes per century = %6ld\n", nb, nc);
      fflush(stdout);
    }
    printf("%16.2f %16.2f %16.2f\n", js->eh_ss[0], js->eh_ss[1], js->eh_ss[2]);
    for (i = 0; i < js->eh_ncon; ++i)
      printf("%.6s\t%24.16f\n", js->ch_cnam + i * 6, js->eh_cval[i]);
    fflush(stdout);
  }
#endif
  return OK;
}

static void reorder(char *x, int size, int number)
{
  int i, j;
  char s[8];
  char *sp1 = x;
  char *sp2 = &s[0];
  for (i = 0; i < number; i++) {
    for (j = 0; j < size; j++)
      *(sp2 + j) = *(sp1 + size - j - 1);
    for (j = 0; j < size; j++)
      *(sp1 + j) = *(sp2 + j);
    sp1 += size;
  }
}

void swi_close_jpl_file(void)
{
  if (js != NULL) {
    if (js->jplfptr != NULL)
      fclose(js->jplfptr);
    if (js->jplfname != NULL)
      FREE((void *) js->jplfname);
    if (js->jplfpath != NULL)
      FREE((void *) js->jplfpath);
    FREE((void *) js);
    js = NULL;
  }
}

int swi_open_jpl_file(double *ss, char *fname, char *fpath, char *serr)
{
  int retc = OK;

  if (js != NULL && js->jplfptr != NULL)
    return OK;
  if ((js = (struct jpl_save *) CALLOC(1, sizeof(struct jpl_save))) == NULL
    || (js->jplfname = (char *) MALLOC(strlen(fname)+1)) == NULL
    || (js->jplfpath = (char *) MALLOC(strlen(fpath)+1)) == NULL
    ) {
    if (serr != NULL)
      strcpy(serr, "error in malloc() with JPL ephemeris.");
    return ERR;
  }
  strcpy(js->jplfname, fname);
  strcpy(js->jplfpath, fpath);
  retc = read_const_jpl(ss, serr);
  if (retc != OK)
    swi_close_jpl_file();
  else {

    js->pc[0] = 1;
    js->pc[1] = 2;
    js->vc[1] = 1;
    js->ac[2] = 4;
    js->jc[3] = 24;
  }
  return retc;
}

int32 swi_get_jpl_denum(void)
{
  return js->eh_denum;
}
