#include "swejpl.h"
#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"
#include <time.h>

#define SEFLG_EPHMASK	(SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)
static int find_maximum(double y00, double y11, double y2, double dx,
			double *dxret, double *yret);
static int find_zero(double y00, double y11, double y2, double dx,
			double *dxret, double *dxret2);
static double calc_dip(double geoalt, double atpress, double attemp, double lapse_rate);
static double calc_astronomical_refr(double geoalt,double atpress, double attemp);
static TLS double const_lapse_rate = SE_LAPSE_RATE;

#if 0
#define DSUN 	(1391978489.9 / AUNIT)

#else
#define DSUN 	(1392000000.0 / AUNIT)
#endif
#define DMOON 	(3476300.0 / AUNIT)
#define DEARTH  (6378140.0 * 2 / AUNIT)
#define RSUN	(DSUN / 2)
#define RMOON	(DMOON / 2)
#define REARTH	(DEARTH / 2)

static int32 eclipse_where( double tjd_ut, int32 ipl, char *starname, int32 ifl, double *geopos,
	double *dcore, char *serr);
static int32 eclipse_how( double tjd_ut, int32 ipl, char *starname, int32 ifl,
        double geolon, double geolat, double geohgt,
	double *attr, char *serr);
static int32 eclipse_when_loc(double tjd_start, int32 ifl, double *geopos,
	double *tret, double *attr, AS_BOOL backward, char *serr);
static int32 occult_when_loc(double tjd_start, int32 ipl, char *starname, int32 ifl,
        double *geopos, double *tret, double *attr, AS_BOOL backward, char *serr);
static int32 lun_eclipse_how(double tjd_ut, int32 ifl, double *attr,
        double *dcore, char *serr);
static int32 calc_mer_trans(
               double tjd_ut, int32 ipl, int32 epheflag, int32 rsmi,
               double *geopos,
               char *starname,
               double *tret,
               char *serr);
static int32 calc_planet_star(double tjd_et, int32 ipl, char *starname, int32 iflag, double *x, char *serr);

struct saros_data {int series_no; double tstart;};

#define SAROS_CYCLE  6585.3213
#define NSAROS_SOLAR 181
struct saros_data saros_data_solar[NSAROS_SOLAR] = {
{0, 641886.5},
{1, 672214.5},
{2, 676200.5},
{3, 693357.5},
{4, 723685.5},
{5, 727671.5},
{6, 744829.5},
{7, 775157.5},
{8, 779143.5},
{9, 783131.5},
{10, 820044.5},
{11, 810859.5},
{12, 748993.5},
{13, 792492.5},
{14, 789892.5},
{15, 787294.5},
{16, 824207.5},
{17, 834779.5},
{18, 838766.5},
{19, 869094.5},
{20, 886251.5},
{21, 890238.5},
{22, 927151.5},
{23, 937722.5},
{24, 941709.5},
{25, 978623.5},
{26, 989194.5},
{27, 993181.5},
{28, 1023510.5},
{29, 1034081.5},
{30, 972214.5},
{31, 1061811.5},
{32, 1006529.5},
{33, 997345.5},
{34, 1021088.5},
{35, 1038245.5},
{36, 1042231.5},
{37, 1065974.5},
{38, 1089716.5},
{39, 1093703.5},
{40, 1117446.5},
{41, 1141188.5},
{42, 1145175.5},
{43, 1168918.5},
{44, 1192660.5},
{45, 1196647.5},
{46, 1220390.5},
{47, 1244132.5},
{48, 1234948.5},
{49, 1265277.5},
{50, 1282433.5},
{51, 1207395.5},
{52, 1217968.5},
{53, 1254881.5},
{54, 1252282.5},
{55, 1262855.5},
{56, 1293182.5},
{57, 1297169.5},
{58, 1314326.5},
{59, 1344654.5},
{60, 1348640.5},
{61, 1365798.5},
{62, 1396126.5},
{63, 1400112.5},
{64, 1417270.5},
{65, 1447598.5},
{66, 1444999.5},
{67, 1462157.5},
{68, 1492485.5},
{69, 1456959.5},
{70, 1421434.5},
{71, 1471518.5},
{72, 1455748.5},
{73, 1466320.5},
{74, 1496648.5},
{75, 1500634.5},
{76, 1511207.5},
{77, 1548120.5},
{78, 1552106.5},
{79, 1562679.5},
{80, 1599592.5},
{81, 1603578.5},
{82, 1614150.5},
{83, 1644479.5},
{84, 1655050.5},
{85, 1659037.5},
{86, 1695950.5},
{87, 1693351.5},
{88, 1631484.5},
{89, 1727666.5},
{90, 1672384.5},
{91, 1663200.5},
{92, 1693529.5},
{93, 1710685.5},
{94, 1714672.5},
{95, 1738415.5},
{96, 1755572.5},
{97, 1766144.5},
{98, 1789887.5},
{99, 1807044.5},
{100, 1817616.5},
{101, 1841359.5},
{102, 1858516.5},
{103, 1862502.5},
{104, 1892831.5},
{105, 1903402.5},
{106, 1887633.5},
{107, 1924547.5},
{108, 1921948.5},
{109, 1873251.5},
{110, 1890409.5},
{111, 1914151.5},
{112, 1918138.5},
{113, 1935296.5},
{114, 1959038.5},
{115, 1963024.5},
{116, 1986767.5},
{117, 2010510.5},
{118, 2014496.5},
{119, 2031654.5},
{120, 2061982.5},
{121, 2065968.5},
{122, 2083126.5},
{123, 2113454.5},
{124, 2104269.5},
{125, 2108256.5},
{126, 2151755.5},
{127, 2083302.5},
{128, 2080704.5},
{129, 2124203.5},
{130, 2121603.5},
{131, 2132176.5},
{132, 2162504.5},
{133, 2166490.5},
{134, 2177062.5},
{135, 2207390.5},
{136, 2217962.5},
{137, 2228534.5},
{138, 2258862.5},
{139, 2269434.5},
{140, 2273421.5},
{141, 2310334.5},
{142, 2314320.5},
{143, 2311722.5},
{144, 2355221.5},
{145, 2319695.5},
{146, 2284169.5},
{147, 2314498.5},
{148, 2325069.5},
{149, 2329056.5},
{150, 2352799.5},
{151, 2369956.5},
{152, 2380528.5},
{153, 2404271.5},
{154, 2421428.5},
{155, 2425414.5},
{156, 2455743.5},
{157, 2472900.5},
{158, 2476886.5},
{159, 2500629.5},
{160, 2517786.5},
{161, 2515187.5},
{162, 2545516.5},
{163, 2556087.5},
{164, 2487635.5},
{165, 2504793.5},
{166, 2535121.5},
{167, 2525936.5},
{168, 2543094.5},
{169, 2573422.5},
{170, 2577408.5},
{171, 2594566.5},
{172, 2624894.5},
{173, 2628880.5},
{174, 2646038.5},
{175, 2669780.5},
{176, 2673766.5},
{177, 2690924.5},
{178, 2721252.5},
{179, 2718653.5},
{180, 2729226.5},
};

#define NSAROS_LUNAR 180
struct saros_data saros_data_lunar[NSAROS_LUNAR] = {
{1, 782437.5},
{2, 799593.5},
{3, 783824.5},
{4, 754884.5},
{5, 824724.5},
{6, 762857.5},
{7, 773430.5},
{8, 810343.5},
{9, 807743.5},
{10, 824901.5},
{11, 855229.5},
{12, 859215.5},
{13, 876373.5},
{14, 906701.5},
{15, 910687.5},
{16, 927845.5},
{17, 958173.5},
{18, 962159.5},
{19, 979317.5},
{20, 1009645.5},
{21, 1007046.5},
{22, 1017618.5},
{23, 1054531.5},
{24, 979493.5},
{25, 976895.5},
{26, 1020394.5},
{27, 1017794.5},
{28, 1028367.5},
{29, 1058695.5},
{30, 1062681.5},
{31, 1073253.5},
{32, 1110167.5},
{33, 1114153.5},
{34, 1131311.5},
{35, 1161639.5},
{36, 1165625.5},
{37, 1176197.5},
{38, 1213111.5},
{39, 1217097.5},
{40, 1221084.5},
{41, 1257997.5},
{42, 1255398.5},
{43, 1186946.5},
{44, 1283128.5},
{45, 1227845.5},
{46, 1225247.5},
{47, 1255575.5},
{48, 1272732.5},
{49, 1276719.5},
{50, 1307047.5},
{51, 1317619.5},
{52, 1328191.5},
{53, 1358519.5},
{54, 1375676.5},
{55, 1379663.5},
{56, 1409991.5},
{57, 1420562.5},
{58, 1424549.5},
{59, 1461463.5},
{60, 1465449.5},
{61, 1436509.5},
{62, 1493179.5},
{63, 1457653.5},
{64, 1435298.5},
{65, 1452456.5},
{66, 1476198.5},
{67, 1480184.5},
{68, 1503928.5},
{69, 1527670.5},
{70, 1531656.5},
{71, 1548814.5},
{72, 1579142.5},
{73, 1583128.5},
{74, 1600286.5},
{75, 1624028.5},
{76, 1628015.5},
{77, 1651758.5},
{78, 1675500.5},
{79, 1672901.5},
{80, 1683474.5},
{81, 1713801.5},
{82, 1645349.5},
{83, 1649336.5},
{84, 1686249.5},
{85, 1683650.5},
{86, 1694222.5},
{87, 1731136.5},
{88, 1735122.5},
{89, 1745694.5},
{90, 1776022.5},
{91, 1786594.5},
{92, 1797166.5},
{93, 1827494.5},
{94, 1838066.5},
{95, 1848638.5},
{96, 1878966.5},
{97, 1882952.5},
{98, 1880354.5},
{99, 1923853.5},
{100, 1881741.5},
{101, 1852801.5},
{102, 1889715.5},
{103, 1893701.5},
{104, 1897688.5},
{105, 1928016.5},
{106, 1938588.5},
{107, 1942575.5},
{108, 1972903.5},
{109, 1990059.5},
{110, 1994046.5},
{111, 2024375.5},
{112, 2034946.5},
{113, 2045518.5},
{114, 2075847.5},
{115, 2086418.5},
{116, 2083820.5},
{117, 2120733.5},
{118, 2124719.5},
{119, 2062852.5},
{120, 2086596.5},
{121, 2103752.5},
{122, 2094568.5},
{123, 2118311.5},
{124, 2142054.5},
{125, 2146040.5},
{126, 2169783.5},
{127, 2186940.5},
{128, 2197512.5},
{129, 2214670.5},
{130, 2238412.5},
{131, 2242398.5},
{132, 2266142.5},
{133, 2289884.5},
{134, 2287285.5},
{135, 2311028.5},
{136, 2334770.5},
{137, 2292659.5},
{138, 2276890.5},
{139, 2326974.5},
{140, 2304619.5},
{141, 2308606.5},
{142, 2345520.5},
{143, 2349506.5},
{144, 2360078.5},
{145, 2390406.5},
{146, 2394392.5},
{147, 2411550.5},
{148, 2441878.5},
{149, 2445864.5},
{150, 2456437.5},
{151, 2486765.5},
{152, 2490751.5},
{153, 2501323.5},
{154, 2538236.5},
{155, 2529052.5},
{156, 2473771.5},
{157, 2563367.5},
{158, 2508085.5},
{159, 2505486.5},
{160, 2542400.5},
{161, 2546386.5},
{162, 2556958.5},
{163, 2587287.5},
{164, 2597858.5},
{165, 2601845.5},
{166, 2632173.5},
{167, 2649330.5},
{168, 2653317.5},
{169, 2683645.5},
{170, 2694217.5},
{171, 2698203.5},
{172, 2728532.5},
{173, 2739103.5},
{174, 2683822.5},
{175, 2740492.5},
{176, 2724722.5},
{177, 2708952.5},
{178, 2732695.5},
{179, 2749852.5},
{180, 2753839.5},
};

int32 CALL_CONV swe_sol_eclipse_where(
		double tjd_ut,
                int32 ifl,
		double *geopos,
		double *attr,
		char *serr)
{
  int32 retflag, retflag2;
  double dcore[10];
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_ut, ifl, 0, serr);
  if ((retflag = eclipse_where(tjd_ut, SE_SUN, NULL, ifl, geopos, dcore, serr)) < 0)
    return retflag;
  if ((retflag2 = eclipse_how(tjd_ut, SE_SUN, NULL, ifl, geopos[0], geopos[1], 0, attr, serr)) == ERR)
    return retflag2;
  attr[3] = dcore[0];
  return retflag;
}

int32 CALL_CONV swe_lun_occult_where(
		double tjd_ut,
                int32 ipl,
                char *starname,
                int32 ifl,
		double *geopos,
		double *attr,
		char *serr)
{
  int32 retflag, retflag2;
  double dcore[10];
  if (ipl < 0) ipl = 0;
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_ut, ifl, 0, serr);

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  if ((retflag = eclipse_where(tjd_ut, ipl, starname, ifl, geopos, dcore, serr)) < 0)
    return retflag;
  if ((retflag2 = eclipse_how(tjd_ut, ipl, starname, ifl, geopos[0], geopos[1], 0, attr, serr)) == ERR)
    return retflag2;
  attr[3] = dcore[0];
  return retflag;
}

static int32 eclipse_where( double tjd_ut, int32 ipl, char *starname, int32 ifl, double *geopos, double *dcore,
		char *serr)
{
  int i;
  int32 retc = 0, niter = 0;
  double e[6], et[6], rm[6], rs[6], rmt[6], rst[6], xs[6], xst[6];
#if 0
  double erm[6];
#endif
  double x[6];
  double lm[6], ls[6], lx[6];
  double dsm, dsmt, d0, D0, s0, r0, d, s, dm;
  double de = 6378140.0 / AUNIT;
  double earthobl = 1 - EARTH_OBLATENESS;
  double deltat, tjd, sidt;
  double drad;
  double sinf1, sinf2, cosf1, cosf2;
  double rmoon = RMOON;
  double dmoon = 2 * rmoon;
  int32 iflag, iflag2;

  AS_BOOL no_eclipse = FALSE;
  struct epsilon *oe = &swed.oec;
  for (i = 0; i < 10; i++)
    dcore[i] = 0;

  iflag = SEFLG_SPEED | SEFLG_EQUATORIAL | ifl;
  iflag2 = iflag | SEFLG_RADIANS;
  iflag  = iflag | SEFLG_XYZ;
  deltat = swe_deltat_ex(tjd_ut, ifl, serr);
  tjd = tjd_ut + deltat;

  if ((retc = swe_calc(tjd, SE_MOON, iflag, rm, serr)) == ERR)
    return retc;

  if ((retc = swe_calc(tjd, SE_MOON, iflag2, lm, serr)) == ERR)
    return retc;

  if ((retc = calc_planet_star(tjd, ipl, starname, iflag, rs, serr)) == ERR)
    return retc;

  if ((retc = calc_planet_star(tjd, ipl, starname, iflag2, ls, serr)) == ERR)
    return retc;

  for (i = 0; i <= 2; i++)
    rst[i] = rs[i];

  for (i = 0; i <= 2; i++)
    rmt[i] = rm[i];
  if (iflag & SEFLG_NONUT)
    sidt = swe_sidtime0(tjd_ut, oe->eps * RADTODEG, 0) * 15 * DEGTORAD;
  else
    sidt = swe_sidtime(tjd_ut) * 15 * DEGTORAD;

  if (starname != NULL && *starname != '\0')
    drad = 0;
  else if (ipl < NDIAM)
    drad = pla_diam[ipl] / 2 / AUNIT;
  else if (ipl > SE_AST_OFFSET)
    drad = swed.ast_diam / 2 * 1000 / AUNIT;
  else
    drad = 0;
iter_where:
  for (i = 0; i <= 2; i++) {
    rs[i] = rst[i];
    rm[i] = rmt[i];
  }

  for (i = 0; i <= 2; i++)
    lx[i] = lm[i];
  swi_polcart(lx, rm);
  rm[2] /= earthobl;

  dm = sqrt(square_sum(rm));

  for (i = 0; i <= 2; i++)
    lx[i] = ls[i];
  swi_polcart(lx, rs);
  rs[2] /= earthobl;

  for (i = 0; i <= 2; i++) {
    e[i] = (rm[i] - rs[i]);
    et[i] = (rmt[i] - rst[i]);
  }

  dsm = sqrt(square_sum(e));
  dsmt = sqrt(square_sum(et));

  for (i = 0; i <= 2; i++) {
    e[i] /= dsm;
    et[i] /= dsmt;
#if 0
    erm[i] = rm[i] / dm;
#endif
  }
  sinf1 = ((drad - rmoon) / dsm);
  cosf1 = sqrt(1 - sinf1 * sinf1);
  sinf2 = ((drad + rmoon) / dsm);
  cosf2 = sqrt(1 - sinf2 * sinf2);

  s0 = -dot_prod(rm, e);

  r0 = sqrt(dm * dm - s0 * s0);

  d0 = (s0 / dsm * (drad * 2 - dmoon) - dmoon) / cosf1;

  D0 = (s0 / dsm * (drad * 2 + dmoon) + dmoon) / cosf2;
  dcore[2] = r0;
  dcore[3] = d0;
  dcore[4] = D0;
  dcore[5] = cosf1;
  dcore[6] = cosf2;
  for (i = 2; i < 5; i++)
    dcore[i] *= AUNIT / 1000.0;

  retc = 0;
  if (de * cosf1 >= r0) {
    retc |= SE_ECL_CENTRAL;
  } else if (r0 <= de * cosf1 + fabs(d0) / 2) {
    retc |= SE_ECL_NONCENTRAL;
  } else if (r0 <= de * cosf2 + D0 / 2) {
    retc |= (SE_ECL_PARTIAL | SE_ECL_NONCENTRAL);
  } else {
    if (serr != NULL)
      sprintf(serr, "no solar eclipse at tjd = %f", tjd);
    for (i = 0; i < 2; i++)
      geopos[i] = 0;
    *dcore = 0;
    retc = 0;
    d = 0;
    no_eclipse = TRUE;

  }

  d = s0 * s0 + de * de - dm * dm;
  if (d > 0)
    d = sqrt(d);
  else
    d = 0;

  s = s0 - d;

#if 0

  if (d == 0) {
    double ds, a, b;

    ds = sqrt(square_sum(rs));
    a = PI - acos(swi_dot_prod_unit(e, erm));

    b = 34.4556 / 60.0 * DEGTORAD + asin(drad / ds);
# if 0

    if (retc & SE_ECL_PARTIAL) {
      d = d0;
      sinf = sinf1;
    } else {
      d = D0;
      sinf = sinf2;
    }
    c = (r0 - de) / d * 2 * sinf;
    if (c > sinf1) {
      b -= .....;
    }
      printf("%f %f %f", a * RADTODEG, b * RADTODEG, s);
      printf(" %f\n", s);
# else
    if (retc & SE_ECL_PARTIAL)
      b -= asin(sinf2);
    else
      b -= asin(sinf1);
# endif
    s += tan(b) * cos(PI / 2 - a) * dm;
  }
#endif

  for (i = 0; i <= 2; i++)
    xs[i] = rm[i] + s * e[i];

  for (i = 0; i <= 2; i++)
    xst[i] = xs[i];
  xst[2] *= earthobl;
  swi_cartpol(xst, xst);
  if (niter <= 0) {
    double cosfi = cos(xst[1]);
    double sinfi = sin(xst[1]);
    double eobl = EARTH_OBLATENESS;
    double cc= 1 / sqrt(cosfi * cosfi + (1-eobl) * (1-eobl) * sinfi * sinfi);
    double ss= (1-eobl) * (1-eobl) * cc;
    earthobl =  ss;
    niter++;
    goto iter_where;
  }
  swi_polcart(xst, xst);

  swi_cartpol(xs, xs);

  xs[0] -= sidt;
  xs[0] *= RADTODEG;
  xs[1] *= RADTODEG;
  xs[0] = swe_degnorm(xs[0]);

  if (xs[0] > 180)
    xs[0] -= 360;
  geopos[0] = xs[0];
  geopos[1] = xs[1];

  for (i = 0; i <= 2; i++)
    x[i] = rmt[i] - xst[i];
  s = sqrt(square_sum(x));

  *dcore = (s / dsmt * ( drad * 2 - dmoon) - dmoon) * cosf1;
  *dcore *= AUNIT / 1000.0;

  dcore[1] = (s / dsmt * ( drad * 2 + dmoon) + dmoon) * cosf2;
  dcore[1] *= AUNIT / 1000.0;
  if (!(retc & SE_ECL_PARTIAL) && !no_eclipse) {
    if (*dcore > 0) {

      retc |= SE_ECL_ANNULAR;
    } else {

      retc |= SE_ECL_TOTAL;
    }
  }
  return retc;
}

static int32 calc_planet_star(double tjd_et, int32 ipl, char *starname, int32 iflag, double *x, char *serr)
{
  int retc = OK;
  if (starname == NULL || *starname == '\0') {
    retc = swe_calc(tjd_et, ipl, iflag, x, serr);
  } else {
    retc = swe_fixstar(starname, tjd_et, iflag, x, serr);
  }
  return retc;
}

int32 CALL_CONV swe_sol_eclipse_how(
          double tjd_ut,
          int32 ifl,
          double *geopos,
          double *attr,
          char *serr)
{
  int32 retflag, retflag2, i;
  double dcore[10], ls[6], xaz[6];
  double geopos2[20];
  for (i = 0; i <= 10; i++)
    attr[i] = 0;
  if (geopos[2] < SEI_ECL_GEOALT_MIN || geopos[2] > SEI_ECL_GEOALT_MAX) {
    if (serr != NULL)
      sprintf(serr, "location for eclipses must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_ut, ifl, 0, serr);
  if ((retflag = eclipse_how(tjd_ut, SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr)) == ERR)
    return retflag;
  if ((retflag2 = eclipse_where(tjd_ut, SE_SUN, NULL, ifl, geopos2, dcore, serr)) == ERR)
    return retflag2;
  if (retflag)
    retflag |= (retflag2 & (SE_ECL_CENTRAL | SE_ECL_NONCENTRAL));
  attr[3] = dcore[0];
  swe_set_topo(geopos[0], geopos[1], geopos[2]);
  if (swe_calc_ut(tjd_ut, SE_SUN, ifl | SEFLG_TOPOCTR | SEFLG_EQUATORIAL, ls, serr) == ERR)
    return ERR;
  swe_azalt(tjd_ut, SE_EQU2HOR, geopos, 0, 10, ls, xaz);
  attr[4] = xaz[0];
  attr[5] = xaz[1];
  attr[6] = xaz[2];
  if (xaz[2] <= 0)
    retflag = 0;
  if (retflag == 0) {
    for (i = 0; i <= 3; i++)
      attr[i] = 0;
    for (i = 8; i <= 10; i++)
      attr[i] = 0;
  }
  return retflag;
}

#define USE_AZ_NAV 0
static int32 eclipse_how( double tjd_ut, int32 ipl, char *starname, int32 ifl,
          double geolon, double geolat, double geohgt,
          double *attr, char *serr)
{
  int i, j, k;
  int32 retc = 0;
  double te, d;
  double xs[6], xm[6], ls[6], lm[6], x1[6], x2[6];
  double rmoon, rsun, rsplusrm, rsminusrm;
  double dctr;
  double drad;
  int32 iflag = SEFLG_EQUATORIAL | SEFLG_TOPOCTR | ifl;
  int32 iflagcart = iflag | SEFLG_XYZ;
#if USE_AZ_NAV
  double mdd, eps, sidt, armc;
#endif
  double xh[6], hmin_appr;
  double lsun, lmoon, lctr, lsunleft, a, b, sc1, sc2;
  double geopos[3];
  for (i = 0; i < 10; i++)
    attr[i] = 0;
  geopos[0] = geolon;
  geopos[1] = geolat;
  geopos[2] = geohgt;
  te = tjd_ut + swe_deltat_ex(tjd_ut, ifl, serr);
  swe_set_topo(geolon, geolat, geohgt);
  if (calc_planet_star(te, ipl, starname, iflag, ls, serr) == ERR)
    return ERR;
  if (swe_calc(te, SE_MOON, iflag, lm, serr) == ERR)
    return ERR;
  if (calc_planet_star(te, ipl, starname, iflagcart, xs, serr) == ERR)
    return ERR;
  if (swe_calc(te, SE_MOON, iflagcart, xm, serr) == ERR)
    return ERR;

  if (starname != NULL && *starname != '\0')
    drad = 0;
  else if (ipl < NDIAM)
    drad = pla_diam[ipl] / 2 / AUNIT;
  else if (ipl > SE_AST_OFFSET)
    drad = swed.ast_diam / 2 * 1000 / AUNIT;
  else
    drad = 0;

#if USE_AZ_NAV
  eps = swi_epsiln(te, iflag);
  if (iflag & SEFLG_NONUT)
    sidt = swe_sidtime0(tjd_ut, eps * RADTODEG, 0) * 15;
  else
    sidt = swe_sidtime(tjd_ut) * 15;
  armc = sidt + geolon;
  mdd = swe_degnorm(ls[0] - armc);
  xh[0] = swe_degnorm(mdd - 90);
  xh[1] = ls[1];
  xh[2] = ls[2];
  swe_cotrans(xh, xh, 90 - geolat);
#else
  swe_azalt(tjd_ut, SE_EQU2HOR, geopos, 0, 10, ls, xh);
#endif

  rmoon = asin(RMOON / lm[2]) * RADTODEG;
  rsun = asin(drad / ls[2]) * RADTODEG;
  rsplusrm = rsun + rmoon;
  rsminusrm = rsun - rmoon;
  for (i = 0; i < 3; i++) {
    x1[i] = xs[i] / ls[2];
    x2[i] = xm[i] / lm[2];
  }
  dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;

  if (dctr < rsminusrm)
    retc = SE_ECL_ANNULAR;
  else if (dctr < fabs(rsminusrm))
    retc = SE_ECL_TOTAL;
  else if (dctr < rsplusrm)
    retc = SE_ECL_PARTIAL;
  else {
    retc = 0;
    if (serr != NULL)
      sprintf(serr, "no solar eclipse at tjd = %f", tjd_ut);
  }

  if (rsun > 0)
    attr[1] = rmoon / rsun;
  else
    attr[1] = 0;

  lsun = asin(rsun / 2 * DEGTORAD) * 2;
  lsunleft = (-dctr + rsun + rmoon);
  if (lsun > 0) {
    attr[0] = lsunleft / rsun / 2;
  } else {

    attr[0] = 1;
  }

  lsun = rsun;
  lmoon = rmoon;
  lctr = dctr;
  if (retc == 0 || lsun == 0) {

    attr[2] = 1;
  } else if (retc == SE_ECL_TOTAL || retc == SE_ECL_ANNULAR) {
    attr[2] = lmoon * lmoon / lsun / lsun;
  } else {
    a = 2 * lctr * lmoon;
    b = 2 * lctr * lsun;
    if (a < 1e-9) {
      attr[2] = lmoon * lmoon / lsun / lsun;
    } else {
      a = (lctr * lctr + lmoon * lmoon - lsun * lsun) / a;
      if (a > 1) a = 1;
      if (a < -1) a = -1;
      b = (lctr * lctr + lsun * lsun - lmoon * lmoon) / b;
      if (b > 1) b = 1;
      if (b < -1) b = -1;
      a = acos(a);
      b = acos(b);
      sc1 = a * lmoon * lmoon / 2;
      sc2 = b * lsun * lsun / 2;
      sc1 -= (cos(a) * sin(a)) * lmoon * lmoon / 2;
      sc2 -= (cos(b) * sin(b)) * lsun * lsun / 2;
      attr[2] = (sc1 + sc2) * 2 / PI / lsun / lsun;
    }
  }
  attr[7] = dctr;

  hmin_appr = -(34.4556 + (1.75 + 0.37) * sqrt(geohgt)) / 60;
  if (xh[1] + rsun + fabs(hmin_appr) >= 0 && retc)
    retc |= SE_ECL_VISIBLE;
#if USE_AZ_NAV
  attr[4] = swe_degnorm(90 - xh[0]);
#else
  attr[4] = xh[0];
#endif
  attr[5] = xh[1];
  attr[6] = xh[2];
  if (ipl == SE_SUN && (starname == NULL || *starname == '\0')) {

    attr[8] = attr[0];
    if (retc & (SE_ECL_TOTAL | SE_ECL_ANNULAR))
      attr[8] = attr[1];

    for (i = 0; i < NSAROS_SOLAR; i++) {
      d = (tjd_ut - saros_data_solar[i].tstart) / SAROS_CYCLE;
      if (d < 0 && d * SAROS_CYCLE > -2) d = 0.0000001;
      if (d < 0) continue;
      j = (int) d;
      if ((d - j) * SAROS_CYCLE < 2) {
	attr[9] = (double) saros_data_solar[i].series_no;
	attr[10] = (double) j + 1;
	break;
      }
      k = j + 1;
      if ((k - d) * SAROS_CYCLE < 2) {
	attr[9] = (double) saros_data_solar[i].series_no;
	attr[10] = (double) k + 1;
	break;
      }
    }
    if (i == NSAROS_SOLAR) {
      attr[9] = attr[10] = -99999999;
    }
  }
  return retc;
}

int32 CALL_CONV swe_sol_eclipse_when_glob(double tjd_start, int32 ifl, int32 ifltype,
     double *tret, int32 backward, char *serr)
{
  int i, j, k, m, n, o, i1 = 0, i2 = 0;
  int32 retflag = 0, retflag2 = 0;
  double de = 6378.140, a;
  double t, tt, tjd, tjds, dt, dtint, dta, dtb;
  double T, T2, T3, T4, K, M, Mm;
  double E, Ff;
  double xs[6], xm[6], ls[6], lm[6];
  double rmoon, rsun, dcore[10];
  double dc[3], dctr;
  double twohr = 2.0 / 24.0;
  double tenmin = 10.0 / 24.0 / 60.0;
  double dt1 = 0, dt2 = 0;
  double geopos[20], attr[20];
  double dtstart, dtdiv;
  double xa[6], xb[6];
  int direction = 1;
  AS_BOOL dont_times = FALSE;
  int32 iflag, iflagcart;
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_start, ifl, 0, serr);
  iflag = SEFLG_EQUATORIAL | ifl;
  iflagcart = iflag | SEFLG_XYZ;
  if (ifltype == (SE_ECL_PARTIAL | SE_ECL_CENTRAL)) {
    if (serr != NULL)
      strcpy(serr, "central partial eclipses do not exist");
    return ERR;
  }
  if (ifltype == (SE_ECL_ANNULAR_TOTAL | SE_ECL_NONCENTRAL)) {
    if (serr != NULL)
      strcpy(serr, "non-central hybrid (annular-total) eclipses do not exist");
    return ERR;
  }
  if (ifltype == 0)
    ifltype = SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_PARTIAL
           | SE_ECL_ANNULAR_TOTAL | SE_ECL_NONCENTRAL | SE_ECL_CENTRAL;
  if (ifltype == SE_ECL_TOTAL || ifltype == SE_ECL_ANNULAR || ifltype == SE_ECL_ANNULAR_TOTAL)
    ifltype |= (SE_ECL_NONCENTRAL | SE_ECL_CENTRAL);
  if (ifltype == SE_ECL_PARTIAL)
    ifltype |= SE_ECL_NONCENTRAL;
  if (backward)
    direction = -1;
  K = (int) ((tjd_start - J2000) / 365.2425 * 12.3685);
  K -= direction;
next_try:
  retflag = 0;
  dont_times = FALSE;
  for (i = 0; i <= 9; i++)
    tret[i] = 0;
  T = K / 1236.85;
  T2 = T * T; T3 = T2 * T; T4 = T3 * T;
  Ff = swe_degnorm(160.7108 + 390.67050274 * K
               - 0.0016341 * T2
               - 0.00000227 * T3
               + 0.000000011 * T4);
  if (Ff > 180)
    Ff -= 180;
  if (Ff > 21 && Ff < 159) {
    K += direction;
    goto next_try;
  }

  tjd = 2451550.09765 + 29.530588853 * K
                      + 0.0001337 * T2
                      - 0.000000150 * T3
                      + 0.00000000073 * T4;
  M = swe_degnorm(2.5534 + 29.10535669 * K
                      - 0.0000218 * T2
                      - 0.00000011 * T3);
  Mm = swe_degnorm(201.5643 + 385.81693528 * K
                      + 0.1017438 * T2
                      + 0.00001239 * T3
                      + 0.000000058 * T4);
  E = 1 - 0.002516 * T - 0.0000074 * T2;
  M *= DEGTORAD;
  Mm *= DEGTORAD;
  tjd = tjd - 0.4075 * sin(Mm)
            + 0.1721 * E * sin(M);

  dtstart = 1;
  if (tjd < 2000000 || tjd > 2500000)
    dtstart = 5;
  dtdiv = 4;
  for (dt = dtstart;
       dt > 0.0001;
       dt /= dtdiv) {
    for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {
      if (swe_calc(t, SE_SUN, iflag, ls, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_SUN, iflagcart, xs, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      for (m = 0; m < 3; m++) {
        xa[m] = xs[m] / ls[2];
        xb[m] = xm[m] / lm[2];
      }
      dc[i] = acos(swi_dot_prod_unit(xa, xb)) * RADTODEG;
      rmoon = asin(RMOON / lm[2]) * RADTODEG;
      rsun = asin(RSUN / ls[2]) * RADTODEG;
      dc[i] -= (rmoon + rsun);
    }
    find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
    tjd += dtint + dt;
  }
  tjds = tjd - swe_deltat_ex(tjd, ifl, serr);
  tjds = tjd - swe_deltat_ex(tjds, ifl, serr);
  tjds = tjd = tjd - swe_deltat_ex(tjds, ifl, serr);
  if ((retflag = eclipse_where(tjd, SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
    return retflag;
  retflag2 = retflag;

  if ((retflag2 = eclipse_how(tjd, SE_SUN, NULL, ifl, geopos[0], geopos[1], 0, attr, serr)) == ERR)
    return retflag2;
  if (retflag2 == 0) {
    K += direction;
    goto next_try;
  }
  tret[0] = tjd;
  if ((backward && tret[0] >= tjd_start - 0.0001)
    || (!backward && tret[0] <= tjd_start + 0.0001)) {
    K += direction;
    goto next_try;
  }

  if ((retflag = eclipse_where(tjd, SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
    return retflag;
  if (retflag == 0) {
    retflag = SE_ECL_PARTIAL | SE_ECL_NONCENTRAL;
    tret[4] = tret[5] = tjd;
    dont_times = TRUE;
  }

  if (!(ifltype & SE_ECL_NONCENTRAL) && (retflag & SE_ECL_NONCENTRAL)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_CENTRAL) && (retflag & SE_ECL_CENTRAL)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_ANNULAR) && (retflag & SE_ECL_ANNULAR)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_PARTIAL) && (retflag & SE_ECL_PARTIAL)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & (SE_ECL_TOTAL | SE_ECL_ANNULAR_TOTAL)) && (retflag & SE_ECL_TOTAL)) {
    K += direction;
    goto next_try;
  }
  if (dont_times)
    goto end_search_global;

  if (retflag & SE_ECL_PARTIAL)
    o = 0;
  else if (retflag & SE_ECL_NONCENTRAL)
    o = 1;
  else
    o = 2;
  dta = twohr;
  dtb = tenmin / 3.0;
  for (n = 0; n <= o; n++) {
    if (n == 0) {

      i1 = 2; i2 = 3;
    } else if (n == 1) {
      if (retflag & SE_ECL_PARTIAL)
        continue;
      i1 = 4; i2 = 5;
    } else if (n == 2) {
      if (retflag & SE_ECL_NONCENTRAL)
        continue;
      i1 = 6; i2 = 7;
    }
    for (i = 0, t = tjd - dta; i <= 2; i += 1, t += dta) {
      if ((retflag2 = eclipse_where(t, SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
        return retflag2;
      if (n == 0)
        dc[i] = dcore[4] / 2 + de / dcore[5] - dcore[2];
      else if (n == 1)
        dc[i] = fabs(dcore[3]) / 2 + de / dcore[6] - dcore[2];
      else if (n == 2)
        dc[i] = de / dcore[6] - dcore[2];
    }
    find_zero(dc[0], dc[1], dc[2], dta, &dt1, &dt2);
    tret[i1] = tjd + dt1 + dta;
    tret[i2] = tjd + dt2 + dta;
    for (m = 0, dt = dtb; m < 3; m++, dt /= 3) {
      for (j = i1; j <= i2; j += (i2 - i1)) {
        for (i = 0, t = tret[j] - dt; i < 2; i++, t += dt) {
          if ((retflag2 = eclipse_where(t, SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
            return retflag2;
          if (n == 0)
            dc[i] = dcore[4] / 2 + de / dcore[5] - dcore[2];
          else if (n == 1)
            dc[i] = fabs(dcore[3]) / 2 + de / dcore[6] - dcore[2];
          else if (n == 2)
            dc[i] = de / dcore[6] - dcore[2];
        }
        dt1 = dc[1] / ((dc[1] - dc[0]) / dt);
        tret[j] -= dt1;
      }
    }
  }

  if (retflag & SE_ECL_TOTAL) {
    if ((retflag2 = eclipse_where(tret[0], SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
      return retflag2;
    dc[0] = *dcore;
    if ((retflag2 = eclipse_where(tret[4], SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
      return retflag2;
    dc[1] = *dcore;
    if ((retflag2 = eclipse_where(tret[5], SE_SUN, NULL, ifl, geopos, dcore, serr)) == ERR)
      return retflag2;
    dc[2] = *dcore;

    if (dc[0] * dc[1] < 0 || dc[0] * dc[2] < 0) {
      retflag |= SE_ECL_ANNULAR_TOTAL;
      retflag &= ~SE_ECL_TOTAL;
    }
  }

  if (!(ifltype & SE_ECL_TOTAL) && (retflag & SE_ECL_TOTAL)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_ANNULAR_TOTAL) && (retflag & SE_ECL_ANNULAR_TOTAL)) {
    K += direction;
    goto next_try;
  }

  k = 2;
  for (i = 0; i < 2; i++) {
    j = i + k;
    tt = tret[j] + swe_deltat_ex(tret[j], ifl, serr);
    if (swe_calc(tt, SE_SUN, iflag, ls, serr) == ERR)
        return ERR;
    if (swe_calc(tt, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
    dc[i] = swe_degnorm(ls[0] - lm[0]);
    if (dc[i] > 180)
      dc[i] -= 360;
  }
  if (dc[0] * dc[1] >= 0)
    tret[1] = 0;
  else {
    tjd = tjds;
    dt = 0.1;
    dt1 = (tret[3] - tret[2]) / 2.0;
    if (dt1 < dt)
      dt = dt1 / 2.0;
    for (j = 0;
        dt > 0.01;
        j++, dt /= 3) {
      for (i = 0, t = tjd; i <= 1; i++, t -= dt) {
        tt = t + swe_deltat_ex(t, ifl, serr);
        if (swe_calc(tt, SE_SUN, iflag, ls, serr) == ERR)
          return ERR;
        if (swe_calc(tt, SE_MOON, iflag, lm, serr) == ERR)
          return ERR;
        dc[i] = swe_degnorm(ls[0] - lm[0]);
        if (dc[i] > 180)
          dc[i] -= 360;
        if (dc[i] > 180)
          dc[i] -= 360;
      }
      a = (dc[1] - dc[0]) / dt;
      if (a < 1e-10)
        break;
      dt1 = dc[0] / a;
      tjd += dt1;
    }
    tret[1] = tjd;
  }
end_search_global:
  return retflag;

}

int32 CALL_CONV swe_lun_occult_when_glob(
     double tjd_start, int32 ipl, char *starname, int32 ifl, int32 ifltype,
     double *tret, int32 backward, char *serr)
{
  int i, j, k, m, n, o, i1 = 0, i2 = 0;
  int32 retflag = 0, retflag2 = 0;
  double de = 6378.140, a;
  double t, tt, tjd = 0, tjds = 0, dt, dtint, dta, dtb;
  double drad, dl;
  double xs[6], xm[6], ls[6], lm[6];
  double rmoon, rsun, dcore[10];
  double dc[20], dctr;
  double twohr = 2.0 / 24.0;
  double tenmin = 10.0 / 24.0 / 60.0;
  double dt1 = 0, dt2 = 0, dadd2 = 1;
  double geopos[20];
  double dtstart, dtdiv;
  int direction = 1;
  int32 ifltype2;
  int32 iflag, iflagcart;
  AS_BOOL dont_times = FALSE;
  int32 one_try = backward & SE_ECL_ONE_TRY;
  if (ipl < 0) ipl = 0;

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_start, ifl, 0, serr);
  iflag = SEFLG_EQUATORIAL | ifl;
  iflagcart = iflag | SEFLG_XYZ;
  backward &= 1L;

  if (ifltype == (SE_ECL_PARTIAL | SE_ECL_CENTRAL)) {
    if (serr != NULL)
      strcpy(serr, "central partial eclipses do not exist");
    return ERR;
  }
  if (ipl != SE_SUN) {
    ifltype2 = (ifltype & ~(SE_ECL_NONCENTRAL | SE_ECL_CENTRAL));
    if (ifltype2 == SE_ECL_ANNULAR || ifltype == SE_ECL_ANNULAR_TOTAL) {
      if (serr != NULL)
	sprintf(serr, "annular occulation do not exist for object %d %s\n", ipl, starname);
      return ERR;
    }
  }
  if (ipl != SE_SUN && (ifltype & (SE_ECL_ANNULAR | SE_ECL_ANNULAR_TOTAL)))
    ifltype &= ~(SE_ECL_ANNULAR|SE_ECL_ANNULAR_TOTAL);
  if (ifltype == 0) {
    ifltype = SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_NONCENTRAL | SE_ECL_CENTRAL;
    if (ipl == SE_SUN)
      ifltype |= (SE_ECL_ANNULAR | SE_ECL_ANNULAR_TOTAL);
  }
  if (ifltype & (SE_ECL_TOTAL | SE_ECL_ANNULAR | SE_ECL_ANNULAR_TOTAL))
    ifltype |= (SE_ECL_NONCENTRAL | SE_ECL_CENTRAL);
  if (ifltype & SE_ECL_PARTIAL)
    ifltype |= SE_ECL_NONCENTRAL;
  retflag = 0;
  for (i = 0; i <= 9; i++)
    tret[i] = 0;
  if (backward)
    direction = -1;
  t = tjd_start;
  tjd = t;
next_try:
  if (calc_planet_star(t, ipl, starname, ifl, ls, serr) == ERR)
      return ERR;

  if (fabs(ls[1]) > 7 && starname != NULL && *starname != '\0') {
    if (serr != NULL)
      sprintf(serr, "occultation never occurs: star %s has ecl. lat. %.1f", starname, ls[1]);
    return ERR;
  }
  if (swe_calc(t, SE_MOON, ifl, lm, serr) == ERR)
      return ERR;
  dl = swe_degnorm(ls[0] - lm[0]);
  if (direction < 0)
    dl -= 360;

  while (fabs(dl) > 0.1) {
    t += dl / 13;
    if (calc_planet_star(t, ipl, starname, ifl, ls, serr) == ERR)
	return ERR;
    if (swe_calc(t, SE_MOON, ifl, lm, serr) == ERR)
	return ERR;
    dl = swe_degnorm(ls[0] - lm[0]);
    if (dl > 180) dl -= 360;
  }
  tjd = t;

  drad = fabs(ls[1] - lm[1]);
  if (drad > 2) {
    if (one_try) {
      tret[0] = t + direction;
      return 0;
    }
    t += direction * 20;
    tjd = t;
    goto next_try;
  }

  if (starname != NULL && *starname != '\0')
    drad = 0;
  else if (ipl < NDIAM)
    drad = pla_diam[ipl] / 2 / AUNIT;
  else if (ipl > SE_AST_OFFSET)
    drad = swed.ast_diam / 2 * 1000 / AUNIT;
  else
    drad = 0;

  dtstart = dadd2;
  dtdiv = 3;
  for (dt = dtstart;
       dt > 0.0001;
       dt /= dtdiv) {
    for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {
      if (calc_planet_star(t, ipl, starname, iflag, ls, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
      if (calc_planet_star(t, ipl, starname, iflagcart, xs, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      dc[i] = acos(swi_dot_prod_unit(xs, xm)) * RADTODEG;
      rmoon = asin(RMOON / lm[2]) * RADTODEG;
      rsun = asin(drad / ls[2]) * RADTODEG;
      dc[i] -= (rmoon + rsun);
    }
    find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
    tjd += dtint + dt;
  }
  tjd -= swe_deltat_ex(tjd, ifl, serr);
  tjds = tjd;
  if ((retflag = eclipse_where(tjd, ipl, starname, ifl, geopos, dcore, serr)) == ERR)
    return retflag;
  retflag2 = retflag;

  if (retflag2 == 0) {

    if (one_try) {
      tret[0] = tjd;
      return 0;
    }

    t = tjd + direction * 20;
    tjd = t;
    goto next_try;
  }
  tret[0] = tjd;

  if ((backward && tret[0] >= tjd_start - 0.0001)
    || (!backward && tret[0] <= tjd_start + 0.0001)) {

    t = tjd + direction * 20;
    tjd = t;
    goto next_try;
  }

  if ((retflag = eclipse_where(tjd, ipl, starname, ifl, geopos, dcore, serr)) == ERR)
    return retflag;
  if (retflag == 0) {
    retflag = SE_ECL_PARTIAL | SE_ECL_NONCENTRAL;
    tret[4] = tret[5] = tjd;
    dont_times = TRUE;
  }

  if (!(ifltype & SE_ECL_NONCENTRAL) && (retflag & SE_ECL_NONCENTRAL)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_CENTRAL) && (retflag & SE_ECL_CENTRAL)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_ANNULAR) && (retflag & SE_ECL_ANNULAR)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_PARTIAL) && (retflag & SE_ECL_PARTIAL)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }

  if (!(ifltype & (SE_ECL_TOTAL | SE_ECL_ANNULAR_TOTAL)) && (retflag & SE_ECL_TOTAL)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }
  if (dont_times)
    goto end_search_global;

  if (retflag & SE_ECL_PARTIAL)
    o = 0;
  else if (retflag & SE_ECL_NONCENTRAL)
    o = 1;
  else
    o = 2;
  dta = twohr;
  dtb = tenmin;
  for (n = 0; n <= o; n++) {
    if (n == 0) {

      i1 = 2; i2 = 3;
    } else if (n == 1) {
      if (retflag & SE_ECL_PARTIAL)
        continue;
      i1 = 4; i2 = 5;
    } else if (n == 2) {
      if (retflag & SE_ECL_NONCENTRAL)
        continue;
      i1 = 6; i2 = 7;
    }
    for (i = 0, t = tjd - dta; i <= 2; i += 1, t += dta) {
      if ((retflag2 = eclipse_where(t, ipl, starname, ifl, geopos, dcore, serr)) == ERR)
        return retflag2;
      if (n == 0)
        dc[i] = dcore[4] / 2 + de / dcore[5] - dcore[2];
      else if (n == 1)
        dc[i] = fabs(dcore[3]) / 2 + de / dcore[6] - dcore[2];
      else if (n == 2)
        dc[i] = de / dcore[6] - dcore[2];
    }
    find_zero(dc[0], dc[1], dc[2], dta, &dt1, &dt2);
    tret[i1] = tjd + dt1 + dta;
    tret[i2] = tjd + dt2 + dta;
    for (m = 0, dt = dtb; m < 3; m++, dt /= 3) {
      for (j = i1; j <= i2; j += (i2 - i1)) {
        for (i = 0, t = tret[j] - dt; i < 2; i++, t += dt) {
          if ((retflag2 = eclipse_where(t, ipl, starname, ifl, geopos, dcore, serr)) == ERR)
            return retflag2;
          if (n == 0)
            dc[i] = dcore[4] / 2 + de / dcore[5] - dcore[2];
          else if (n == 1)
            dc[i] = fabs(dcore[3]) / 2 + de / dcore[6] - dcore[2];
          else if (n == 2)
            dc[i] = de / dcore[6] - dcore[2];
        }
        dt1 = dc[1] / ((dc[1] - dc[0]) / dt);
        tret[j] -= dt1;
      }
    }
  }

  if (retflag & SE_ECL_TOTAL) {
    if ((retflag2 = eclipse_where(tret[0], ipl, starname, ifl, geopos, dcore, serr)) == ERR)
      return retflag2;
    dc[0] = *dcore;
    if ((retflag2 = eclipse_where(tret[4], ipl, starname, ifl, geopos, dcore, serr)) == ERR)
      return retflag2;
    dc[1] = *dcore;
    if ((retflag2 = eclipse_where(tret[5], ipl, starname, ifl, geopos, dcore, serr)) == ERR)
      return retflag2;
    dc[2] = *dcore;

    if (dc[0] * dc[1] < 0 || dc[0] * dc[2] < 0) {
      retflag |= SE_ECL_ANNULAR_TOTAL;
      retflag &= ~SE_ECL_TOTAL;
    }
  }

  if (!(ifltype & SE_ECL_TOTAL) && (retflag & SE_ECL_TOTAL)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_ANNULAR_TOTAL) && (retflag & SE_ECL_ANNULAR_TOTAL)) {

    t = tjd + direction * 20;
    if (one_try) {
      tret[0] = tjd;
      return 0;
    }
    tjd = t;
    goto next_try;
  }

  k = 2;
  for (i = 0; i < 2; i++) {
    j = i + k;
    tt = tret[j] + swe_deltat_ex(tret[j], ifl, serr);
    if (calc_planet_star(tt, ipl, starname, iflag, ls, serr) == ERR)
        return ERR;
    if (swe_calc(tt, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
    dc[i] = swe_degnorm(ls[0] - lm[0]);
    if (dc[i] > 180)
      dc[i] -= 360;
  }
  if (dc[0] * dc[1] >= 0)
    tret[1] = 0;
  else {
    tjd = tjds;
    dt = 0.1;
    dt1 = (tret[3] - tret[2]) / 2.0;
    if (dt1 < dt)
      dt = dt1 / 2.0;
    for (j = 0;
        dt > 0.01;
        j++, dt /= 3) {
      for (i = 0, t = tjd; i <= 1; i++, t -= dt) {
        tt = t + swe_deltat_ex(t, ifl, serr);
        if (calc_planet_star(tt, ipl, starname, iflag, ls, serr) == ERR)
          return ERR;
        if (swe_calc(tt, SE_MOON, iflag, lm, serr) == ERR)
          return ERR;
        dc[i] = swe_degnorm(ls[0] - lm[0]);
        if (dc[i] > 180)
          dc[i] -= 360;
        if (dc[i] > 180)
          dc[i] -= 360;
      }
      a = (dc[1] - dc[0]) / dt;
      if (a < 1e-10)
        break;
      dt1 = dc[0] / a;
      tjd += dt1;
    }
    tret[1] = tjd;
  }
end_search_global:
  return retflag;

}

int32 CALL_CONV swe_sol_eclipse_when_loc(double tjd_start, int32 ifl,
     double *geopos, double *tret, double *attr, int32 backward, char *serr)
{
  int32 retflag = 0, retflag2 = 0;
  double geopos2[20], dcore[10];
  if (geopos[2] < SEI_ECL_GEOALT_MIN || geopos[2] > SEI_ECL_GEOALT_MAX) {
    if (serr != NULL)
      sprintf(serr, "location for eclipses must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_start, ifl, 0, serr);
  if ((retflag = eclipse_when_loc(tjd_start, ifl, geopos, tret, attr, backward, serr)) <= 0)
    return retflag;

  if ((retflag2 = eclipse_where(tret[0], SE_SUN, NULL, ifl, geopos2, dcore, serr)) == ERR)
    return retflag2;
  retflag |= (retflag2 & SE_ECL_NONCENTRAL);
  attr[3] = dcore[0];
  return retflag;
}

int32 CALL_CONV swe_lun_occult_when_loc(double tjd_start, int32 ipl, char *starname, int32 ifl,
     double *geopos, double *tret, double *attr, int32 backward, char *serr)
{
  int32 retflag = 0, retflag2 = 0;
  double geopos2[20], dcore[10];

  if (geopos[2] < SEI_ECL_GEOALT_MIN || geopos[2] > SEI_ECL_GEOALT_MAX) {
    if (serr != NULL)
      sprintf(serr, "location for occultations must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  if (ipl < 0) ipl = 0;
  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_start, ifl, 0, serr);
  if ((retflag = occult_when_loc(tjd_start, ipl, starname, ifl, geopos, tret, attr, backward, serr)) <= 0)
    return retflag;

  if ((retflag2 = eclipse_where(tret[0], ipl, starname, ifl, geopos2, dcore, serr)) == ERR)
    return retflag2;
  retflag |= (retflag2 & SE_ECL_NONCENTRAL);
  attr[3] = dcore[0];
  return retflag;
}

static int32 eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr)
{
  int i, j, k, m;
  int32 retflag = 0, retc;
  double t, tjd, dt, dtint, K, T, T2, T3, T4, F, M, Mm;
  double tjdr, tjds = 0;
  double E, Ff;
  double xs[6], xm[6], ls[6], lm[6], x1[6], x2[6], dm, ds;
  double rmoon, rsun, rsplusrm, rsminusrm;
  double dc[3], dctr, dctrmin;
  double twomin = 2.0 / 24.0 / 60.0;
  double tensec = 10.0 / 24.0 / 60.0 / 60.0;
  double twohr = 2.0 / 24.0;
  double tenmin = 10.0 / 24.0 / 60.0;
  double dt1 = 0, dt2 = 0, dtdiv, dtstart;
  int32 iflag = SEFLG_EQUATORIAL | SEFLG_TOPOCTR | ifl;
  int32 iflagcart = iflag | SEFLG_XYZ;
  swe_set_topo(geopos[0], geopos[1], geopos[2]);
  K = (int) ((tjd_start - J2000) / 365.2425 * 12.3685);
  if (backward)
    K++;
  else
    K--;
next_try:
  T = K / 1236.85;
  T2 = T * T; T3 = T2 * T; T4 = T3 * T;
  Ff = F = swe_degnorm(160.7108 + 390.67050274 * K
               - 0.0016341 * T2
               - 0.00000227 * T3
               + 0.000000011 * T4);
  if (Ff > 180)
    Ff -= 180;
  if (Ff > 21 && Ff < 159) {
    if (backward)
      K--;
    else
      K++;
    goto next_try;
  }

  tjd = 2451550.09765 + 29.530588853 * K
                      + 0.0001337 * T2
                      - 0.000000150 * T3
                      + 0.00000000073 * T4;
  M = swe_degnorm(2.5534 + 29.10535669 * K
                      - 0.0000218 * T2
                      - 0.00000011 * T3);
  Mm = swe_degnorm(201.5643 + 385.81693528 * K
                      + 0.1017438 * T2
                      + 0.00001239 * T3
                      + 0.000000058 * T4);

  E = 1 - 0.002516 * T - 0.0000074 * T2;

  M *= DEGTORAD;
  Mm *= DEGTORAD;
  F *= DEGTORAD;

  tjd = tjd - 0.4075 * sin(Mm)
            + 0.1721 * E * sin(M);
  swe_set_topo(geopos[0], geopos[1], geopos[2]);
  dtdiv = 2;
  dtstart = 0.5;
  if (tjd < 1900000 || tjd > 2500000)
    dtstart = 2;
  for (dt = dtstart;
       dt > 0.00001;
       dt /= dtdiv) {
    if (dt < 0.1)
      dtdiv = 3;
    for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {

      if (swe_calc(t, SE_SUN, iflagcart, xs, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_SUN, iflag, ls, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
      dm = sqrt(square_sum(xm));
      ds = sqrt(square_sum(xs));
      for (k = 0; k < 3; k++) {
        x1[k] = xs[k] / ds ;
        x2[k] = xm[k] / dm ;
      }
      dc[i] = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
    }
    find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
    tjd += dtint + dt;
  }
  if (swe_calc(tjd, SE_SUN, iflagcart, xs, serr) == ERR)
        return ERR;
  if (swe_calc(tjd, SE_SUN, iflag, ls, serr) == ERR)
        return ERR;
  if (swe_calc(tjd, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
  if (swe_calc(tjd, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
  dctr = acos(swi_dot_prod_unit(xs, xm)) * RADTODEG;
  rmoon = asin(RMOON / lm[2]) * RADTODEG;
  rsun = asin(RSUN / ls[2]) * RADTODEG;
  rsplusrm = rsun + rmoon;
  rsminusrm = rsun - rmoon;
  if (dctr > rsplusrm) {
    if (backward)
      K--;
    else
      K++;
    goto next_try;
  }
  tret[0] = tjd - swe_deltat_ex(tjd, ifl, serr);
  tret[0] = tjd - swe_deltat_ex(tret[0], ifl, serr);
  if ((backward && tret[0] >= tjd_start - 0.0001)
    || (!backward && tret[0] <= tjd_start + 0.0001)) {
    if (backward)
      K--;
    else
      K++;
    goto next_try;
  }
  if (dctr < rsminusrm)
    retflag = SE_ECL_ANNULAR;
  else if (dctr < fabs(rsminusrm))
    retflag = SE_ECL_TOTAL;
  else if (dctr <= rsplusrm)
    retflag = SE_ECL_PARTIAL;
  dctrmin = dctr;

  if (dctr > fabs(rsminusrm))
    tret[2] = tret[3] = 0;
  else {
    dc[1] = fabs(rsminusrm) - dctrmin;
    for (i = 0, t = tjd - twomin; i <= 2; i += 2, t = tjd + twomin) {
      if (swe_calc(t, SE_SUN, iflagcart, xs, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      dm = sqrt(square_sum(xm));
      ds = sqrt(square_sum(xs));
      rmoon = asin(RMOON / dm) * RADTODEG;
      rmoon *= 0.99916;
      rsun = asin(RSUN / ds) * RADTODEG;
      rsminusrm = rsun - rmoon;
      for (k = 0; k < 3; k++) {
        x1[k] = xs[k] / ds ;
        x2[k] = xm[k] / dm ;
      }
      dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
      dc[i] = fabs(rsminusrm) - dctr;
    }
    find_zero(dc[0], dc[1], dc[2], twomin, &dt1, &dt2);
    tret[2] = tjd + dt1 + twomin;
    tret[3] = tjd + dt2 + twomin;
    for (m = 0, dt = tensec; m < 2; m++, dt /= 10) {
      for (j = 2; j <= 3; j++) {
        if (swe_calc(tret[j], SE_SUN, iflagcart | SEFLG_SPEED, xs, serr) == ERR)
          return ERR;
        if (swe_calc(tret[j], SE_MOON, iflagcart | SEFLG_SPEED, xm, serr) == ERR)
          return ERR;
        for (i = 0; i < 2; i++) {
          if (i == 1) {
            for(k = 0; k < 3; k++) {
              xs[k] -= xs[k+3] * dt;
              xm[k] -= xm[k+3] * dt;
            }
          }
          dm = sqrt(square_sum(xm));
          ds = sqrt(square_sum(xs));
          rmoon = asin(RMOON / dm) * RADTODEG;
	  rmoon *= 0.99916;
          rsun = asin(RSUN / ds) * RADTODEG;
          rsminusrm = rsun - rmoon;
          for (k = 0; k < 3; k++) {
            x1[k] = xs[k] / ds ;
            x2[k] = xm[k] / dm ;
          }
          dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
          dc[i] = fabs(rsminusrm) - dctr;
        }
        dt1 = -dc[0] / ((dc[0] - dc[1]) / dt);
        tret[j] += dt1;
      }
    }
    tret[2] -= swe_deltat_ex(tret[2], ifl, serr);
    tret[3] -= swe_deltat_ex(tret[3], ifl, serr);
  }

  dc[1] = rsplusrm - dctrmin;
  for (i = 0, t = tjd - twohr; i <= 2; i += 2, t = tjd + twohr) {
    if (swe_calc(t, SE_SUN, iflagcart, xs, serr) == ERR)
      return ERR;
    if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
      return ERR;
    dm = sqrt(square_sum(xm));
    ds = sqrt(square_sum(xs));
    rmoon = asin(RMOON / dm) * RADTODEG;
    rsun = asin(RSUN / ds) * RADTODEG;
    rsplusrm = rsun + rmoon;
    for (k = 0; k < 3; k++) {
      x1[k] = xs[k] / ds ;
      x2[k] = xm[k] / dm ;
    }
    dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
    dc[i] = rsplusrm - dctr;
  }
  find_zero(dc[0], dc[1], dc[2], twohr, &dt1, &dt2);
  tret[1] = tjd + dt1 + twohr;
  tret[4] = tjd + dt2 + twohr;
  for (m = 0, dt = tenmin; m < 3; m++, dt /= 10) {
    for (j = 1; j <= 4; j += 3) {
      if (swe_calc(tret[j], SE_SUN, iflagcart | SEFLG_SPEED, xs, serr) == ERR)
        return ERR;
      if (swe_calc(tret[j], SE_MOON, iflagcart | SEFLG_SPEED, xm, serr) == ERR)
        return ERR;
      for (i = 0; i < 2; i++) {
        if (i == 1) {
          for(k = 0; k < 3; k++) {
            xs[k] -= xs[k+3] * dt;
            xm[k] -= xm[k+3] * dt;
          }
        }
        dm = sqrt(square_sum(xm));
        ds = sqrt(square_sum(xs));
        rmoon = asin(RMOON / dm) * RADTODEG;
        rsun = asin(RSUN / ds) * RADTODEG;
        rsplusrm = rsun + rmoon;
        for (k = 0; k < 3; k++) {
          x1[k] = xs[k] / ds ;
          x2[k] = xm[k] / dm ;
        }
        dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
        dc[i] = fabs(rsplusrm) - dctr;
      }
      dt1 = -dc[0] / ((dc[0] - dc[1]) / dt);
      tret[j] += dt1;
    }
  }
  tret[1] -= swe_deltat_ex(tret[1], ifl, serr);
  tret[4] -= swe_deltat_ex(tret[4], ifl, serr);

  for (i = 4; i >= 0; i--) {
    if (tret[i] == 0)
      continue;
    if (eclipse_how(tret[i], SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2],
		attr, serr) == ERR)
      return ERR;

    if (attr[6] > 0) {
      retflag |= SE_ECL_VISIBLE;
      switch(i) {
      case 0: retflag |= SE_ECL_MAX_VISIBLE; break;
      case 1: retflag |= SE_ECL_1ST_VISIBLE; break;
      case 2: retflag |= SE_ECL_2ND_VISIBLE; break;
      case 3: retflag |= SE_ECL_3RD_VISIBLE; break;
      case 4: retflag |= SE_ECL_4TH_VISIBLE; break;
      default:  break;
      }
    }
  }
#if 1
  if (!(retflag & SE_ECL_VISIBLE)) {
    if (backward)
      K--;
    else
      K++;
    goto next_try;
  }
#endif
  if ((retc = swe_rise_trans(tret[1] - 0.001, SE_SUN, NULL, iflag, SE_CALC_RISE|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjdr, serr)) == ERR)
    return ERR;
  if (retc == -2)
    return retflag;
  if ((retc = swe_rise_trans(tret[1] - 0.001, SE_SUN, NULL, iflag, SE_CALC_SET|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjds, serr)) == ERR)
    return ERR;
  if (retc == -2)
    return retflag;
  if (tjds < tret[1] || (tjds > tjdr && tjdr > tret[4])) {
    if (backward)
      K--;
    else
      K++;
    goto next_try;
  }
  if (tjdr > tret[1] && tjdr < tret[4]) {
    tret[5] = tjdr;
    if (!(retflag & SE_ECL_MAX_VISIBLE)) {
      tret[0] = tjdr;
      if ((retc = eclipse_how(tret[5], SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr)) == ERR)
	return ERR;
      retflag &= ~(SE_ECL_TOTAL|SE_ECL_ANNULAR|SE_ECL_PARTIAL);
      retflag |= (retc & (SE_ECL_TOTAL|SE_ECL_ANNULAR|SE_ECL_PARTIAL));
    }
  }
  if (tjds > tret[1] && tjds < tret[4]) {
    tret[6] = tjds;
    if (!(retflag & SE_ECL_MAX_VISIBLE)) {
      tret[0] = tjds;
      if ((retc = eclipse_how(tret[6], SE_SUN, NULL, ifl, geopos[0], geopos[1], geopos[2], attr, serr)) == ERR)
	return ERR;
      retflag &= ~(SE_ECL_TOTAL|SE_ECL_ANNULAR|SE_ECL_PARTIAL);
      retflag |= (retc & (SE_ECL_TOTAL|SE_ECL_ANNULAR|SE_ECL_PARTIAL));
    }
  }
  return retflag;
}

static int32 occult_when_loc(
     double tjd_start, int32 ipl, char *starname,
     int32 ifl, double *geopos, double *tret, double *attr,
     int32 backward, char *serr)
{
  int i, j, k, m;
  int32 retflag = 0, retc;
  double t, tjd, dt, dtint;
  double tjdr, tjds = 0;
  double xs[6], xm[6], ls[6], lm[6], x1[6], x2[6], dm, ds;
  double rmoon, rsun, rsplusrm, rsminusrm;
  double dc[20], dctr, dctrmin;
  double twomin = 2.0 / 24.0 / 60.0;
  double tensec = 10.0 / 24.0 / 60.0 / 60.0;
  double twohr = 2.0 / 24.0;
  double tenmin = 10.0 / 24.0 / 60.0;
  double dt1 = 0, dt2 = 0, dtdiv, dtstart;
  double dadd2 = 1;
  double drad, dl;

  int32 iflag = SEFLG_TOPOCTR | ifl;
  int32 iflaggeo = iflag & ~SEFLG_TOPOCTR;
  int32 iflagcart = iflag | SEFLG_XYZ;
  int direction = 1;
  int32 one_try = backward & SE_ECL_ONE_TRY;
  AS_BOOL stop_after_this = FALSE;
  backward &= 1L;
  retflag = 0;
  swe_set_topo(geopos[0], geopos[1], geopos[2]);
  for (i = 0; i <= 9; i++)
    tret[i] = 0;
  if (backward)
    direction = -1;
  t = tjd_start;
  tjd = tjd_start;
next_try:

  if (calc_planet_star(t, ipl, starname, iflaggeo, ls, serr) == ERR)
      return ERR;

  if (fabs(ls[1]) > 7 && starname != NULL && *starname != '\0') {
    if (serr != NULL)
      sprintf(serr, "occultation never occurs: star %s has ecl. lat. %.1f", starname, ls[1]);
    return ERR;
  }
  if (swe_calc(t, SE_MOON, iflaggeo, lm, serr) == ERR)
      return ERR;
  dl = swe_degnorm(ls[0] - lm[0]);
  if (direction < 0)
    dl -= 360;

  while (fabs(dl) > 0.1) {
    t += dl / 13;
    if (calc_planet_star(t, ipl, starname, iflaggeo, ls, serr) == ERR)
	return ERR;
    if (swe_calc(t, SE_MOON, iflaggeo, lm, serr) == ERR)
	return ERR;
    dl = swe_degnorm(ls[0] - lm[0]);
    if (dl > 180) dl -= 360;
  }
  tjd = t;

  drad = fabs(ls[1] - lm[1]);
  if (drad > 2) {
    if (one_try) {
      tret[0] = t + direction;
      return 0;
    }
    t += direction * 20;
    tjd = t;
    goto next_try;
  }

  if (starname != NULL && *starname != '\0')
    drad = 0;
  else if (ipl < NDIAM)
    drad = pla_diam[ipl] / 2 / AUNIT;
  else if (ipl > SE_AST_OFFSET)
    drad = swed.ast_diam / 2 * 1000 / AUNIT;
  else
    drad = 0;

  dtdiv = 2;
  dtstart = dadd2;
  for (dt = dtstart;
       dt > 0.00001;
       dt /= dtdiv) {
    if (dt < 0.01)
      dtdiv = 2;
    for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {

      if (calc_planet_star(t, ipl, starname, iflagcart, xs, serr) == ERR)
        return ERR;
      if (calc_planet_star(t, ipl, starname, iflag, ls, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
      if (dt < 0.1 && fabs(ls[1] - lm[1]) > 2) {
        if (one_try || stop_after_this) {
          stop_after_this = TRUE;
        } else {

          t = tjd + direction * 20;
          tjd = t;
          goto next_try;
        }
      }
      dc[i] = acos(swi_dot_prod_unit(xs, xm)) * RADTODEG;
      rmoon = asin(RMOON / lm[2]) * RADTODEG;
      rsun = asin(drad / ls[2]) * RADTODEG;
      dc[i] -= (rmoon + rsun);
    }
    find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
    tjd += dtint + dt;
  }
  if (stop_after_this) {
    tret[0] = tjd + direction;
    return 0;
  }
  if (calc_planet_star(tjd, ipl, starname, iflagcart, xs, serr) == ERR)
        return ERR;
  if (calc_planet_star(tjd, ipl, starname, iflag, ls, serr) == ERR)
        return ERR;
  if (swe_calc(tjd, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
  if (swe_calc(tjd, SE_MOON, iflag, lm, serr) == ERR)
        return ERR;
  dctr = acos(swi_dot_prod_unit(xs, xm)) * RADTODEG;
  rmoon = asin(RMOON / lm[2]) * RADTODEG;
  rsun = asin(drad / ls[2]) * RADTODEG;
  rsplusrm = rsun + rmoon;
  rsminusrm = rsun - rmoon;
  if (dctr > rsplusrm) {
    if (one_try) {
      tret[0] = tjd + direction;
      return 0;
    }

    t = tjd + direction * 20;
    tjd = t;
    goto next_try;
  }
  tret[0] = tjd - swe_deltat_ex(tjd, ifl, serr);
  tret[0] = tjd - swe_deltat_ex(tret[0], ifl, serr);
  if ((backward && tret[0] >= tjd_start - 0.0001)
    || (!backward && tret[0] <= tjd_start + 0.0001)) {

    if (one_try) {
      tret[0] = tjd + direction;
      return 0;
    }
    t = tjd + direction * 20;
    tjd = t;
    goto next_try;
  }
  if (dctr < rsminusrm)
    retflag = SE_ECL_ANNULAR;
  else if (dctr < fabs(rsminusrm))
    retflag = SE_ECL_TOTAL;
  else if (dctr <= rsplusrm)
    retflag = SE_ECL_PARTIAL;
  dctrmin = dctr;

  if (dctr > fabs(rsminusrm)) {
    tret[2] = tret[3] = 0;

  } else {
    dc[1] = fabs(rsminusrm) - dctrmin;
    for (i = 0, t = tjd - twomin; i <= 2; i += 2, t = tjd + twomin) {
      if (calc_planet_star(t, ipl, starname, iflagcart, xs, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      dm = sqrt(square_sum(xm));
      ds = sqrt(square_sum(xs));
      rmoon = asin(RMOON / dm) * RADTODEG;
      rmoon *= 0.99916;
      rsun = asin(drad / ds) * RADTODEG;
      rsminusrm = rsun - rmoon;
      for (k = 0; k < 3; k++) {
        x1[k] = xs[k] / ds ;
        x2[k] = xm[k] / dm ;
      }
      dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
      dc[i] = fabs(rsminusrm) - dctr;
    }
    find_zero(dc[0], dc[1], dc[2], twomin, &dt1, &dt2);
    tret[2] = tjd + dt1 + twomin;
    tret[3] = tjd + dt2 + twomin;
    for (m = 0, dt = tensec; m < 2; m++, dt /= 10) {
      for (j = 2; j <= 3; j++) {
        if (calc_planet_star(tret[j], ipl, starname, iflagcart | SEFLG_SPEED, xs, serr) == ERR)
          return ERR;
        if (swe_calc(tret[j], SE_MOON, iflagcart | SEFLG_SPEED, xm, serr) == ERR)
          return ERR;
        for (i = 0; i < 2; i++) {
          if (i == 1) {
            for(k = 0; k < 3; k++) {
              xs[k] -= xs[k+3] * dt;
              xm[k] -= xm[k+3] * dt;
            }
          }
          dm = sqrt(square_sum(xm));
          ds = sqrt(square_sum(xs));
          rmoon = asin(RMOON / dm) * RADTODEG;
	  rmoon *= 0.99916;
          rsun = asin(drad / ds) * RADTODEG;
          rsminusrm = rsun - rmoon;
          for (k = 0; k < 3; k++) {
            x1[k] = xs[k] / ds ;
            x2[k] = xm[k] / dm ;
          }
          dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
          dc[i] = fabs(rsminusrm) - dctr;
        }
        dt1 = -dc[0] / ((dc[0] - dc[1]) / dt);
        tret[j] += dt1;
      }
    }
    tret[2] -= swe_deltat_ex(tret[2], ifl, serr);
    tret[3] -= swe_deltat_ex(tret[3], ifl, serr);

  }

  dc[1] = rsplusrm - dctrmin;
if (starname == NULL || *starname == '\0') {
  for (i = 0, t = tjd - twohr; i <= 2; i += 2, t = tjd + twohr) {
    if (calc_planet_star(t, ipl, starname, iflagcart, xs, serr) == ERR)
      return ERR;
    if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
      return ERR;
    dm = sqrt(square_sum(xm));
    ds = sqrt(square_sum(xs));
    rmoon = asin(RMOON / dm) * RADTODEG;
    rsun = asin(drad / ds) * RADTODEG;
    rsplusrm = rsun + rmoon;
    for (k = 0; k < 3; k++) {
      x1[k] = xs[k] / ds ;
      x2[k] = xm[k] / dm ;
    }
    dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
    dc[i] = rsplusrm - dctr;
  }
  find_zero(dc[0], dc[1], dc[2], twohr, &dt1, &dt2);
  tret[1] = tjd + dt1 + twohr;
  tret[4] = tjd + dt2 + twohr;
  for (m = 0, dt = tenmin; m < 3; m++, dt /= 10) {
    for (j = 1; j <= 4; j += 3) {
      if (calc_planet_star(tret[j], ipl, starname, iflagcart | SEFLG_SPEED, xs, serr) == ERR)
        return ERR;
      if (swe_calc(tret[j], SE_MOON, iflagcart | SEFLG_SPEED, xm, serr) == ERR)
        return ERR;
      for (i = 0; i < 2; i++) {
        if (i == 1) {
          for(k = 0; k < 3; k++) {
            xs[k] -= xs[k+3] * dt;
            xm[k] -= xm[k+3] * dt;
          }
        }
        dm = sqrt(square_sum(xm));
        ds = sqrt(square_sum(xs));
        rmoon = asin(RMOON / dm) * RADTODEG;
        rsun = asin(drad / ds) * RADTODEG;
        rsplusrm = rsun + rmoon;
        for (k = 0; k < 3; k++) {
          x1[k] = xs[k] / ds ;
          x2[k] = xm[k] / dm ;
        }
        dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;
        dc[i] = fabs(rsplusrm) - dctr;
      }
      dt1 = -dc[0] / ((dc[0] - dc[1]) / dt);
      tret[j] += dt1;
    }
  }
  tret[1] -= swe_deltat_ex(tret[1], ifl, serr);
  tret[4] -= swe_deltat_ex(tret[4], ifl, serr);
} else {
  tret[1] = tret[2];
  tret[4] = tret[3];
}

  for (i = 4; i >= 0; i--) {
    if (tret[i] == 0)
      continue;
    if (eclipse_how(tret[i], ipl, starname, ifl, geopos[0], geopos[1], geopos[2],
		attr, serr) == ERR)
      return ERR;

    if (attr[6] > 0) {
      retflag |= SE_ECL_VISIBLE;
      switch(i) {
      case 0: retflag |= SE_ECL_MAX_VISIBLE; break;
      case 1: retflag |= SE_ECL_1ST_VISIBLE; break;
      case 2: retflag |= SE_ECL_2ND_VISIBLE; break;
      case 3: retflag |= SE_ECL_3RD_VISIBLE; break;
      case 4: retflag |= SE_ECL_4TH_VISIBLE; break;
      default:  break;
      }
    }
  }
#if 1
  if (!(retflag & SE_ECL_VISIBLE)) {

    if (one_try) {
      tret[0] = tjd + direction;
      return 0;
    }
    t = tjd + direction * 20;
    tjd = t;
    goto next_try;
  }
#endif
  if ((retc = swe_rise_trans(tret[1] - 0.1, ipl, starname, iflag, SE_CALC_RISE|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjdr, serr)) == ERR)
    return ERR;
  if (retc >= 0 && (retc = swe_rise_trans(tret[1] - 0.1, ipl, starname, iflag, SE_CALC_SET|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjds, serr)) == ERR)
    return ERR;
  if (retc >= 0) {
    if (tjdr > tret[1] && tjdr < tret[4])
      tret[5] = tjdr;
    if (tjds > tret[1] && tjds < tret[4])
      tret[6] = tjds;
  }

    if ((retc = swe_rise_trans(tret[1], SE_SUN, NULL, iflag, SE_CALC_RISE, geopos, 0, 0, &tjdr, serr)) == ERR)
      return ERR;
    if (retc >= 0 && (retc = swe_rise_trans(tret[1], SE_SUN, NULL, iflag, SE_CALC_SET, geopos, 0, 0, &tjds, serr)) == ERR)
      return ERR;
    if (retc >= 0) {
      if (tjds < tjdr)
	retflag |= SE_ECL_OCC_BEG_DAYLIGHT;
    }
    if ((retc = swe_rise_trans(tret[4], SE_SUN, NULL, iflag, SE_CALC_RISE, geopos, 0, 0, &tjdr, serr)) == ERR)
      return ERR;
    if (retc >= 0 && (retc = swe_rise_trans(tret[4], SE_SUN, NULL, iflag, SE_CALC_SET, geopos, 0, 0, &tjds, serr)) == ERR)
      return ERR;
    if (retc >= 0) {
      if (tjds < tjdr)
	retflag |= SE_ECL_OCC_END_DAYLIGHT;
    }

  return retflag;
}

void CALL_CONV swe_azalt(
      double tjd_ut,
      int32  calc_flag,
      double *geopos,
      double atpress,
      double attemp,
      double *xin,
      double *xaz)
{
  int i;
  double x[6], xra[3];
  double armc = swe_degnorm(swe_sidtime(tjd_ut) * 15 + geopos[0]);
  double mdd, eps_true;
  for (i = 0; i < 2; i++)
    xra[i] = xin[i];
  xra[2] = 1;
  if (calc_flag == SE_ECL2HOR) {
    swe_calc(tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL), SE_ECL_NUT, 0, x, NULL);
    eps_true = x[0];
    swe_cotrans(xra, xra, -eps_true);
  }
  mdd = swe_degnorm(xra[0] - armc);
  x[0] = swe_degnorm(mdd - 90);
  x[1] = xra[1];
  x[2] = 1;

  swe_cotrans(x, x, 90 - geopos[1]);

  x[0] = swe_degnorm(x[0] + 90);
  xaz[0] = 360 - x[0];
  xaz[1] = x[1];
  if (atpress == 0) {

    atpress = 1013.25 * pow(1 - 0.0065 * geopos[2] / 288, 5.255);
  }
  xaz[2] = swe_refrac_extended(x[1], geopos[2], atpress, attemp, const_lapse_rate, SE_TRUE_TO_APP, NULL);

}

void CALL_CONV swe_azalt_rev(
      double tjd_ut,
      int32  calc_flag,
      double *geopos,
      double *xin,
      double *xout)
{
  int i;
  double x[6], xaz[3];
  double geolon = geopos[0];
  double geolat = geopos[1];
  double armc = swe_degnorm(swe_sidtime(tjd_ut) * 15 + geolon);
  double eps_true, dang;
  for (i = 0; i < 2; i++)
    xaz[i] = xin[i];
  xaz[2] = 1;

  xaz[0] = 360 - xaz[0];
  xaz[0] = swe_degnorm(xaz[0] - 90);

  dang = geolat - 90;
  swe_cotrans(xaz, xaz, dang);
  xaz[0] = swe_degnorm(xaz[0] + armc + 90);
  xout[0] = xaz[0];
  xout[1] = xaz[1];

  if (calc_flag == SE_HOR2ECL) {
    swe_calc(tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL), SE_ECL_NUT, 0, x, NULL);
    eps_true = x[0];
	swe_cotrans(xaz, x, eps_true);
    xout[0] = x[0];
    xout[1] = x[1];
  }
}

double CALL_CONV swe_refrac(double inalt, double atpress, double attemp, int32 calc_flag)
{
  double a, refr;
  double pt_factor = atpress / 1010.0 * 283.0 / (273.0 + attemp);
  double trualt, appalt;
#if 0

  double y, yy0, D0, N, D, P, Q;
  int i;
  if (calc_flag == SE_TRUE_TO_APP) {
    trualt = inalt;
    if( (trualt < -2.0) || (trualt >= 90.0) )
      return(trualt);

    if( trualt > 15.0 ) {
	  D = 0.00452*atpress/((273.0+attemp)*tan( DEGTORAD*trualt ));
	  return(trualt + D);
	}

    y = trualt;
    D = 0.0;

    P = (atpress - 80.0)/930.0;
    Q = 4.8e-3 * (attemp - 10.0);
    yy0 = y;
    D0 = D;
    for( i=0; i<4; i++ ) {
  	  N = y + (7.31/(y+4.4));
	  N = 1.0/tan(DEGTORAD*N);
	  D = N*P/(60.0 + Q * (N + 39.0));
	  N = y - yy0;
	  yy0 = D - D0 - N;
	  if( (N != 0.0) && (yy0 != 0.0) )

	    N = y - N*(trualt + D - y)/yy0;
	  else

	    N = trualt + D;
	  yy0 = y;
	  D0 = D;
	  y = N;
	}
    return( trualt + D );
  } else {
#else

  if (calc_flag == SE_TRUE_TO_APP) {
    trualt = inalt;
    if (trualt > 15) {
      a = tan((90 - trualt) * DEGTORAD);
      refr = (58.276 * a - 0.0824 * a * a * a);
      refr *=  pt_factor / 3600.0;
    } else if (trualt > -5) {

      a = trualt + 10.3 / (trualt + 5.11);
      if (a + 1e-10 >= 90)
    	refr = 0;
      else
        refr = 1.02 / tan(a * DEGTORAD);
      refr *= pt_factor / 60.0;
    } else
      refr = 0;
    appalt = trualt;
    if (appalt + refr > 0)
      appalt += refr;
    return appalt;
  } else {
#endif

    appalt = inalt;

    a = appalt + 7.31 / (appalt + 4.4);
    if (a + 1e-10 >= 90)
      refr = 0;
    else {
      refr = 1.00 / tan(a * DEGTORAD);
      refr -= 0.06 * sin(14.7 * refr + 13);
    }
    refr *= pt_factor / 60.0;
    trualt = appalt;
    if (appalt - refr > 0)
      trualt = appalt - refr;
    return trualt;
  }
}

void CALL_CONV swe_set_lapse_rate(double lapse_rate)
{
  const_lapse_rate = lapse_rate;
}

double CALL_CONV swe_refrac_extended(double inalt, double geoalt, double atpress, double attemp, double lapse_rate, int32 calc_flag, double *dret)
{
  double refr;
  double trualt;
  double dip = calc_dip(geoalt, atpress, attemp, lapse_rate);
  double D, D0, N, y, yy0;
  int i;

  if( (inalt>90) )
    inalt=180-inalt;
  if (calc_flag == SE_TRUE_TO_APP) {
    if (inalt < -10) {
      if (dret != NULL) {
	dret[0]=inalt;
	dret[1]=inalt;
	dret[2]=0;
	dret[3]=dip;
      }
      return inalt;
    }

    y = inalt;
    D = 0.0;
    yy0 = 0;
    D0 = D;
    for(i=0; i<5; i++) {
      D = calc_astronomical_refr(y,atpress,attemp);
      N = y - yy0;
      yy0 = D - D0 - N;
      if (N != 0.0 && yy0 != 0.0)
        N = y - N*(inalt + D - y)/yy0;
      else
        N = inalt + D;
      yy0 = y;
      D0 = D;
      y = N;
    }
    refr = D;
    if (inalt + refr < dip) {
      if (dret != NULL) {
	dret[0]=inalt;
	dret[1]=inalt;
	dret[2]=0;
	dret[3]=dip;
      }
      return inalt;
    }
    if (dret != NULL) {
      dret[0]=inalt;
      dret[1]=inalt+refr;
      dret[2]=refr;
      dret[3]=dip;
    }
    return inalt+refr;
  } else {
    refr = calc_astronomical_refr(inalt,atpress,attemp);
    trualt=inalt-refr;

    if (dret != NULL) {
      if (inalt > dip) {
	dret[0]=trualt;
	dret[1]=inalt;
	dret[2]=refr;
	dret[3]=dip;
      } else {
	dret[0]=inalt;
	dret[1]=inalt;
	dret[2]=0;
	dret[3]=dip;
      }
    }

    if (inalt >= dip)
      return trualt;
    else
      return inalt;
  }
}

static double calc_astronomical_refr(double inalt,double atpress, double attemp)
{
#if 0

  double refractaccent = 1/tan(DEGTORAD*(inalt + 7.31/(inalt+4.4)));
  double r = (refractaccent - 0.06 * sin(DEGTORAD*(14.7*refractaccent +13)));
  r = ((atpress - 80) / 930 / (1 + 0.00008 * (r + 39) * (attemp - 10)) * r)/60;
 return r;
#else

  double r;
  if (inalt > 17.904104638432) {
    r = 0.97 / tan(inalt * DEGTORAD);
  } else {
    r = (34.46 + 4.23 * inalt + 0.004 * inalt * inalt) / (1 + 0.505 * inalt + 0.0845 * inalt * inalt);
  }
  r = ((atpress - 80) / 930 / (1 + 0.00008 * (r + 39) * (attemp - 10)) * r) / 60.0;
  return r;
#endif
}

static double calc_dip(double geoalt, double atpress, double attemp, double lapse_rate)
{

  double krefr = (0.0342 + lapse_rate) / (0.154 * 0.0238);
  double d = 1-1.8480*krefr*atpress/(273.15+attemp)/(273.15+attemp);

  return -180.0/PI * acos(1 / (1 + geoalt / EARTH_RADIUS)) * sqrt(d);
}

int32 CALL_CONV swe_lun_eclipse_how(
          double tjd_ut,
          int32 ifl,
          double *geopos,
          double *attr,
          char *serr)
{
  double dcore[10];
  double lm[6], xaz[6];
  int32 retc;

  if (geopos != NULL)
    geopos[0] = geopos[0];
  if (geopos != NULL && (geopos[2] < SEI_ECL_GEOALT_MIN || geopos[2] > SEI_ECL_GEOALT_MAX)) {
    if (serr != NULL)
      sprintf(serr, "location for eclipses must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  ifl = ifl & ~SEFLG_TOPOCTR;
  ifl &= ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);
  swi_set_tid_acc(tjd_ut, ifl, 0, serr);
  retc = lun_eclipse_how(tjd_ut, ifl, attr, dcore, serr);
  if (geopos == NULL) {
    return retc;
  }

  swe_set_topo(geopos[0], geopos[1], geopos[2]);
  if (swe_calc_ut(tjd_ut, SE_MOON, ifl | SEFLG_TOPOCTR | SEFLG_EQUATORIAL, lm, serr) == ERR)
    return ERR;
  swe_azalt(tjd_ut, SE_EQU2HOR, geopos, 0, 10, lm, xaz);
  attr[4] = xaz[0];
  attr[5] = xaz[1];
  attr[6] = xaz[2];
  if (xaz[2] <= 0)
    retc = 0;
  return retc;
}

static int32 lun_eclipse_how(
          double tjd_ut,
          int32 ifl,
          double *attr,
          double *dcore,
          char *serr)
{
  int i, j, k;
  int32 retc = 0;
  double e[6], rm[6], rs[6];
  double dsm, d0, D0, s0, r0, ds, dm;
  double dctr, x1[6], x2[6];
  double f1, f2;
  double deltat, tjd, d;
  double cosf1, cosf2;
  double rmoon = RMOON;
  double dmoon = 2 * rmoon;
  int32 iflag;
  for (i = 0; i < 10; i++)
    dcore[i] = 0;
  for (i = 0; i < 20; i++)
    attr[i] = 0;

  iflag = SEFLG_SPEED | SEFLG_EQUATORIAL | ifl;
  iflag  = iflag | SEFLG_XYZ;
  deltat = swe_deltat_ex(tjd_ut, ifl, serr);
  tjd = tjd_ut + deltat;

  if (swe_calc(tjd, SE_MOON, iflag, rm, serr) == ERR)
    return ERR;

  dm = sqrt(square_sum(rm));

  if (swe_calc(tjd, SE_SUN, iflag, rs, serr) == ERR)
    return ERR;

  ds = sqrt(square_sum(rs));
  for (i = 0; i < 3; i++) {
    x1[i] = rs[i] / ds;
    x2[i] = rm[i] / dm;
  }
  dctr = acos(swi_dot_prod_unit(x1, x2)) * RADTODEG;

  for (i = 0; i <= 2; i++)
    rs[i] -= rm[i];

  for (i = 0; i <= 2; i++)
    rm[i] = -rm[i];

  for (i = 0; i <= 2; i++)
    e[i] = (rm[i] - rs[i]);

  dsm = sqrt(square_sum(e));

  for (i = 0; i <= 2; i++)
    e[i] /= dsm;
  f1 = ((RSUN - REARTH) / dsm);
  cosf1 = sqrt(1 - f1 * f1);
  f2 = ((RSUN + REARTH) / dsm);
  cosf2 = sqrt(1 - f2 * f2);

  s0 = -dot_prod(rm, e);

  r0 = sqrt(dm * dm - s0 * s0);

  d0 = fabs(s0 / dsm * (DSUN - DEARTH) - DEARTH) * (1 + 1.0 / 50.0) / cosf1;

  D0 = (s0 / dsm * (DSUN + DEARTH) + DEARTH) * (1 + 1.0 / 50.0) / cosf2;
  d0 /= cosf1;
  D0 /= cosf2;

  d0 *= 0.99405;
  D0 *= 0.98813;
  dcore[0] = r0;
  dcore[1] = d0;
  dcore[2] = D0;
  dcore[3] = cosf1;
  dcore[4] = cosf2;

  retc = 0;
  if (d0 / 2 >= r0 + rmoon / cosf1) {
    retc = SE_ECL_TOTAL;
    attr[0] = (d0 / 2 - r0 + rmoon) / dmoon;
  } else if (d0 / 2 >= r0 - rmoon / cosf1) {
    retc = SE_ECL_PARTIAL;
    attr[0] = (d0 / 2 - r0 + rmoon) / dmoon;
  } else if (D0 / 2 >= r0 - rmoon / cosf2) {
    retc = SE_ECL_PENUMBRAL;
    attr[0] = 0;
  } else {
    if (serr != NULL)
      sprintf(serr, "no lunar eclipse at tjd = %f", tjd);
  }
  attr[8] = attr[0];

  attr[1] = (D0 / 2 - r0 + rmoon) / dmoon;
  if (retc != 0)
    attr[7] = 180 - fabs(dctr);

  for (i = 0; i < NSAROS_LUNAR; i++) {
    d = (tjd_ut - saros_data_lunar[i].tstart) / SAROS_CYCLE;
    if (d < 0 && d * SAROS_CYCLE > -2) d = 0.0000001;
    if (d < 0) continue;
    j = (int) d;
    if ((d - j) * SAROS_CYCLE < 2) {
      attr[9] = (double) saros_data_lunar[i].series_no;
      attr[10] = (double) j + 1;
      break;
    }
    k = j + 1;
    if ((k - d) * SAROS_CYCLE < 2) {
      attr[9] = (double) saros_data_lunar[i].series_no;
      attr[10] = (double) k + 1;
      break;
    }
  }
  if (i == NSAROS_LUNAR) {
    attr[9] = attr[10] = -99999999;
  }
  return retc;
}

int32 CALL_CONV swe_lun_eclipse_when(double tjd_start, int32 ifl, int32 ifltype,
     double *tret, int32 backward, char *serr)
{
  int i, j, m, n, o, i1 = 0, i2 = 0;
  int32 retflag = 0, retflag2 = 0;
  double t, tjd, tjd2, dt, dtint, dta, dtb;
  double T, T2, T3, T4, K, F, M, Mm;
  double E, Ff, F1, A1, Om;
  double xs[6], xm[6], dm, ds;
  double rsun, rearth, dcore[10];
  double dc[3], dctr;
  double twohr = 2.0 / 24.0;
  double tenmin = 10.0 / 24.0 / 60.0;
  double dt1 = 0, dt2 = 0;
  double kk;
  double attr[20];
  double dtstart, dtdiv;
  double xa[6], xb[6];
  int direction = 1;
  int32 iflag;
  int32 iflagcart;
  ifl &= SEFLG_EPHMASK;
  swi_set_tid_acc(tjd_start, ifl, 0, serr);
  iflag = SEFLG_EQUATORIAL | ifl;
  iflagcart = iflag | SEFLG_XYZ;
  ifltype &= ~(SE_ECL_CENTRAL|SE_ECL_NONCENTRAL);
  if (ifltype & (SE_ECL_ANNULAR|SE_ECL_ANNULAR_TOTAL)) {
    ifltype &= ~(SE_ECL_ANNULAR|SE_ECL_ANNULAR_TOTAL);
    if (ifltype == 0) {
      if (serr != NULL) {
        strcpy(serr, "annular lunar eclipses don't exist");
      }
      return ERR;
    }
  }
  if (ifltype == 0)
    ifltype = SE_ECL_TOTAL | SE_ECL_PENUMBRAL | SE_ECL_PARTIAL;
  if (backward)
    direction = -1;
  K = (int) ((tjd_start - J2000) / 365.2425 * 12.3685);
  K -= direction;
next_try:
  retflag = 0;
  for (i = 0; i <= 9; i++)
    tret[i] = 0;
  kk = K + 0.5;
  T = kk / 1236.85;
  T2 = T * T; T3 = T2 * T; T4 = T3 * T;
  Ff = F = swe_degnorm(160.7108 + 390.67050274 * kk
               - 0.0016341 * T2
               - 0.00000227 * T3
               + 0.000000011 * T4);
  if (Ff > 180)
    Ff -= 180;
  if (Ff > 21 && Ff < 159) {
    K += direction;
    goto next_try;
  }

  tjd = 2451550.09765 + 29.530588853 * kk
                      + 0.0001337 * T2
                      - 0.000000150 * T3
                      + 0.00000000073 * T4;
  M = swe_degnorm(2.5534 + 29.10535669 * kk
                      - 0.0000218 * T2
                      - 0.00000011 * T3);
  Mm = swe_degnorm(201.5643 + 385.81693528 * kk
                      + 0.1017438 * T2
                      + 0.00001239 * T3
                      + 0.000000058 * T4);
  Om = swe_degnorm(124.7746 - 1.56375580 * kk
                      + 0.0020691 * T2
                      + 0.00000215 * T3);
  E = 1 - 0.002516 * T - 0.0000074 * T2;
  A1 = swe_degnorm(299.77 + 0.107408 * kk - 0.009173 * T2);
  M *= DEGTORAD;
  Mm *= DEGTORAD;
  F *= DEGTORAD;
  Om *= DEGTORAD;
  F1 = F - 0.02665 * sin(Om) * DEGTORAD;
  A1 *= DEGTORAD;
  tjd = tjd - 0.4075 * sin(Mm)
            + 0.1721 * E * sin(M)
            + 0.0161 * sin(2 * Mm)
            - 0.0097 * sin(2 * F1)
            + 0.0073 * E * sin(Mm - M)
            - 0.0050 * E * sin(Mm + M)
            - 0.0023 * sin(Mm - 2 * F1)
            + 0.0021 * E * sin(2 * M)
            + 0.0012 * sin(Mm + 2 * F1)
            + 0.0006 * E * sin(2 * Mm + M)
            - 0.0004 * sin(3 * Mm)
            - 0.0003 * E * sin(M + 2 * F1)
            + 0.0003 * sin(A1)
            - 0.0002 * E * sin(M - 2 * F1)
            - 0.0002 * E * sin(2 * Mm - M)
            - 0.0002 * sin(Om);

  dtstart = 0.1;
  if (tjd < 2100000 || tjd > 2500000)
    dtstart = 5;
  dtdiv = 4;
  for (j = 0, dt = dtstart;
       dt > 0.001;
       j++, dt /= dtdiv) {
    for (i = 0, t = tjd - dt; i <= 2; i++, t += dt) {
      if (swe_calc(t, SE_SUN, iflagcart, xs, serr) == ERR)
        return ERR;
      if (swe_calc(t, SE_MOON, iflagcart, xm, serr) == ERR)
        return ERR;
      for (m = 0; m < 3; m++) {
        xs[m] -= xm[m];
        xm[m] = -xm[m];
      }
      ds = sqrt(square_sum(xs));
      dm = sqrt(square_sum(xm));
      for (m = 0; m < 3; m++) {
        xa[m] = xs[m] / ds;
        xb[m] = xm[m] / dm;
      }
      dc[i] = acos(swi_dot_prod_unit(xa, xb)) * RADTODEG;
      rearth = asin(REARTH / dm) * RADTODEG;
      rsun = asin(RSUN / ds) * RADTODEG;
      dc[i] -= (rearth + rsun);
    }
    find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dctr);
    tjd += dtint + dt;
  }
  tjd2 = tjd - swe_deltat_ex(tjd, ifl, serr);
  tjd2 = tjd - swe_deltat_ex(tjd2, ifl, serr);
  tjd = tjd - swe_deltat_ex(tjd2, ifl, serr);
  if ((retflag = swe_lun_eclipse_how(tjd, ifl, NULL, attr, serr)) == ERR)
    return retflag;
  if (retflag == 0) {
    K += direction;
    goto next_try;
  }
  tret[0] = tjd;
  if ((backward && tret[0] >= tjd_start - 0.0001)
    || (!backward && tret[0] <= tjd_start + 0.0001)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_PENUMBRAL) && (retflag & SE_ECL_PENUMBRAL)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & SE_ECL_PARTIAL) && (retflag & SE_ECL_PARTIAL)) {
    K += direction;
    goto next_try;
  }

  if (!(ifltype & (SE_ECL_TOTAL)) && (retflag & SE_ECL_TOTAL)) {
    K += direction;
    goto next_try;
  }

  if (retflag & SE_ECL_PENUMBRAL)
    o = 0;
  else if (retflag & SE_ECL_PARTIAL)
    o = 1;
  else
    o = 2;
  dta = twohr;
  dtb = tenmin;
  for (n = 0; n <= o; n++) {
    if (n == 0) {
      i1 = 6; i2 = 7;
    } else if (n == 1) {
      i1 = 2; i2 = 3;
    } else if (n == 2) {
      i1 = 4; i2 = 5;
    }
#if 1
    for (i = 0, t = tjd - dta; i <= 2; i += 1, t += dta) {
      if ((retflag2 = lun_eclipse_how(t, ifl, attr, dcore, serr)) == ERR)
        return retflag2;
      if (n == 0)
        dc[i] = dcore[2] / 2 + RMOON / dcore[4] - dcore[0];
      else if (n == 1)
        dc[i] = dcore[1] / 2 + RMOON / dcore[3] - dcore[0];
      else if (n == 2)
        dc[i] = dcore[1] / 2 - RMOON / dcore[3] - dcore[0];
    }
    find_zero(dc[0], dc[1], dc[2], dta, &dt1, &dt2);
    dtb = (dt1 + dta) / 2;
    tret[i1] = tjd + dt1 + dta;
    tret[i2] = tjd + dt2 + dta;
#else
    tret[i1] = tjd - dtb;
    tret[i2] = tjd + dtb;
#endif
    for (m = 0, dt = dtb / 2; m < 3; m++, dt /= 2) {
      for (j = i1; j <= i2; j += (i2 - i1)) {
        for (i = 0, t = tret[j] - dt; i < 2; i++, t += dt) {
          if ((retflag2 = lun_eclipse_how(t, ifl, attr, dcore, serr)) == ERR)
            return retflag2;
          if (n == 0)
            dc[i] = dcore[2] / 2 + RMOON / dcore[4] - dcore[0];
          else if (n == 1)
            dc[i] = dcore[1] / 2 + RMOON / dcore[3] - dcore[0];
          else if (n == 2)
            dc[i] = dcore[1] / 2 - RMOON / dcore[3] - dcore[0];
        }
        dt1 = dc[1] / ((dc[1] - dc[0]) / dt);
        tret[j] -= dt1;
      }
    }
  }
  return retflag;
}

int32 CALL_CONV swe_lun_eclipse_when_loc(double tjd_start, int32 ifl,
     double *geopos, double *tret, double *attr, int32 backward, char *serr)
{
  int32 retflag = 0, retflag2 = 0, retc;
  double tjdr, tjds = 0, tjd_max = 0;
  int i;
  if (geopos != NULL && (geopos[2] < SEI_ECL_GEOALT_MIN || geopos[2] > SEI_ECL_GEOALT_MAX)) {
    if (serr != NULL)
      sprintf(serr, "location for eclipses must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  ifl &= ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);
next_lun_ecl:
  if ((retflag = swe_lun_eclipse_when(tjd_start, ifl, 0, tret, backward, serr)) == ERR) {
    return ERR;
  }

  retflag = 0;
  for (i = 7; i >= 0; i--) {
    if (i == 1) continue;
    if (tret[i] == 0) continue;
    if ((retflag2 = swe_lun_eclipse_how(tret[i], ifl, geopos, attr, serr)) == ERR)
      return ERR;
    if (attr[6] > 0) {
      retflag |= SE_ECL_VISIBLE;
      switch(i) {
      case 0: retflag |= SE_ECL_MAX_VISIBLE; break;
      case 2: retflag |= SE_ECL_PARTBEG_VISIBLE; break;
      case 3: retflag |= SE_ECL_PARTEND_VISIBLE; break;
      case 4: retflag |= SE_ECL_TOTBEG_VISIBLE; break;
      case 5: retflag |= SE_ECL_TOTEND_VISIBLE; break;
      case 6: retflag |= SE_ECL_PENUMBBEG_VISIBLE; break;
      case 7: retflag |= SE_ECL_PENUMBEND_VISIBLE; break;
      default:  break;
      }
    }
  }
  if (!(retflag & SE_ECL_VISIBLE)) {
    if (backward)
      tjd_start = tret[0] - 25;
    else
      tjd_start = tret[0] + 25;
    goto next_lun_ecl;
  }

  tjd_max = tret[0];
  if ((retc = swe_rise_trans(tret[6] - 0.001, SE_MOON, NULL, ifl, SE_CALC_RISE|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjdr, serr)) == ERR)
    return ERR;
  if (retc >= 0 && (retc = swe_rise_trans(tret[6] - 0.001, SE_MOON, NULL, ifl, SE_CALC_SET|SE_BIT_DISC_BOTTOM, geopos, 0, 0, &tjds, serr)) == ERR)
    return ERR;
  if (retc >= 0) {
    if (tjds < tret[6] || (tjds > tjdr && tjdr > tret[7])) {
      if (backward)
	tjd_start = tret[0] - 25;
      else
	tjd_start = tret[0] + 25;
      goto next_lun_ecl;
    }
    if (tjdr > tret[6] && tjdr < tret[7]) {
      tret[6] = 0;
      for (i = 2; i <= 5; i++) {
	if (tjdr > tret[i])
	  tret[i] = 0;
      }
      tret[8] = tjdr;
      if (tjdr > tret[0]) {
	tjd_max = tjdr;
      }
    }
    if (tjds > tret[6] && tjds < tret[7]) {
      tret[7] = 0;
      for (i = 2; i <= 5; i++) {
	if (tjds < tret[i])
	  tret[i] = 0;
      }
      tret[9] = tjds;
      if (tjds < tret[0]) {
	tjd_max = tjds;
      }
    }
  }
  tret[0] = tjd_max;
  if ((retflag2 = swe_lun_eclipse_how(tjd_max, ifl, geopos, attr, serr)) == ERR)
    return ERR;
  if (retflag2 == 0) {
    if (backward)
      tjd_start = tret[0] - 25;
    else
      tjd_start = tret[0] + 25;
    goto next_lun_ecl;
  }
  retflag |= (retflag2 & SE_ECL_ALLTYPES_LUNAR);
  return retflag;
}

#define EULER 2.718281828459
#define NMAG_ELEM  (SE_VESTA + 1)
#define MAG_MALLAMA_2018  1
#define MAG_MOON_VREIJS   1

static const double mag_elem[NMAG_ELEM][4] = {

                {-26.86, 0, 0, 0},
                {-12.55, 0, 0, 0},

                {-0.42, 3.80, -2.73, 2.00},
                {-4.40, 0.09, 2.39, -0.65},
                {- 1.52, 1.60, 0, 0},
                {- 9.40, 0.5, 0, 0},
                {- 8.88, -2.60, 1.25, 0.044},
                {- 7.19, 0.0, 0, 0},
                {- 6.87, 0.0, 0, 0},
                {- 1.00, 0.0, 0, 0},
                {99, 0, 0, 0},
                {99, 0, 0, 0},
                {99, 0, 0, 0},
                {99, 0, 0, 0},
                {99, 0, 0, 0},

                {6.5, 0.15, 0, 0},
                {7.0, 0.15, 0, 0},
                {3.34, 0.12, 0, 0},
                {4.13, 0.11, 0, 0},
                {5.33, 0.32, 0, 0},
                {3.20, 0.32, 0, 0},
                };
int32 CALL_CONV swe_pheno(double tjd, int32 ipl, int32 iflag, double *attr, char *serr)
{
  int i;
  double xx[6], xx2[6], xxs[6], lbr[6], lbr2[6], dt = 0, dd;
  double fac;
  double T, in, om, sinB;
  double ph1, ph2, me[2];
  int32 iflagp, epheflag, retflag, epheflag2;
  char serr2[AS_MAXCH];
  *serr2 = '\0';
  iflag &= ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  for (i = 0; i < 20; i++)
    attr[i] = 0;

  if (ipl > SE_AST_OFFSET && ipl <= SE_AST_OFFSET + 4)
        ipl = ipl - SE_AST_OFFSET - 1 + SE_CERES;
  iflag = iflag & (SEFLG_EPHMASK |
                   SEFLG_TRUEPOS |
                   SEFLG_J2000 |
                   SEFLG_NONUT |
                   SEFLG_NOGDEFL |
                   SEFLG_NOABERR |
                   SEFLG_TOPOCTR);
  iflagp = iflag & (SEFLG_EPHMASK |
                   SEFLG_TRUEPOS |
                   SEFLG_J2000 |
                   SEFLG_NONUT |
                   SEFLG_NOABERR);
  iflagp |= SEFLG_HELCTR;
  epheflag = iflag & SEFLG_EPHMASK;

  if ((retflag = swe_calc(tjd, (int) ipl, iflag | SEFLG_XYZ, xx, serr)) == ERR)

    return ERR;

  epheflag2 = retflag & SEFLG_EPHMASK;
  if (epheflag != epheflag2) {
    iflag &= ~epheflag;
    iflagp &= ~epheflag;
    iflag |= epheflag2;
    iflagp |= epheflag2;
    epheflag = epheflag2;
  }
  if (swe_calc(tjd, (int) ipl, iflag, lbr, serr) == ERR)

    return ERR;

  if (ipl == SE_MOON) {
    if (swe_calc(tjd, SE_SUN, iflag | SEFLG_XYZ, xxs, serr) == ERR)
      return ERR;
  }
  if (ipl != SE_SUN && ipl != SE_EARTH &&
    ipl != SE_MEAN_NODE && ipl != SE_TRUE_NODE &&
    ipl != SE_MEAN_APOG && ipl != SE_OSCU_APOG) {

    dt = lbr[2] * AUNIT / CLIGHT / 86400.0;
        if (iflag & SEFLG_TRUEPOS)
      dt = 0;

    if (swe_calc(tjd - dt, (int) ipl, iflagp | SEFLG_XYZ, xx2, serr) == ERR)

      return ERR;
    if (swe_calc(tjd - dt, (int) ipl, iflagp, lbr2, serr) == ERR)

      return ERR;

    attr[0] = acos(swi_dot_prod_unit(xx, xx2)) * RADTODEG;

    attr[1] = (1 + cos(attr[0] * DEGTORAD)) / 2;
  }

  if (ipl < NDIAM)
    dd = pla_diam[ipl];
  else if (ipl > SE_AST_OFFSET)
    dd = swed.ast_diam * 1000;
  else
    dd = 0;
  if (lbr[2] < dd / 2 / AUNIT)
    attr[3] = 180;
  else
    attr[3] = asin(dd / 2 / AUNIT / lbr[2]) * 2 * RADTODEG;

  if (ipl > SE_AST_OFFSET || (ipl < NMAG_ELEM && mag_elem[ipl][0] < 99)) {
    if (ipl == SE_SUN) {

      fac = attr[3] / (asin(pla_diam[SE_SUN] / 2.0 / AUNIT) * 2 * RADTODEG);
      fac *= fac;
      attr[4] = mag_elem[ipl][0] - 2.5 * log10(fac);
    } else if (ipl == SE_MOON) {
#if MAG_MOON_VREIJS

     double a = attr[0];
     if (a<=147.1385465) {

       attr[4] = -21.62 + 0.026 * fabs(a) + 0.000000004 * pow(a, 4);
       attr[4]+=5 * log10(lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS);
     } else {

       attr[4] = -4.5444 - (2.5 * log10(pow(180 - a, 3)));
       attr[4]+=5 * log10(lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS);
     }
#else

      attr[4] = -21.62 + 5 * log10(lbr[2] * lbr2[2] * AUNIT / EARTH_RADIUS) + 0.026 * fabs(attr[0]) + 0.000000004 * pow(attr[0], 4);
#endif
#if MAG_MALLAMA_2018

    } else if (ipl == SE_MERCURY) {
      double a = attr[0];
      double a2 = a * a; double a3 = a2 * a; double a4 = a3 * a; double a5 = a4 * a; double a6 = a5 * a;
      attr[4] = -0.613 + a * 6.3280E-02 - a2 * 1.6336E-03 + a3 * 3.3644E-05 - a4 * 3.4265E-07 + a5 * 1.6893E-09 - a6 * 3.0334E-12;
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);
    } else if (ipl == SE_VENUS) {
      double a = attr[0];
      double a2 = a * a; double a3 = a2 * a; double a4 = a3 * a;
      if (a <= 163.7)
	attr[4] = -4.384 - a * 1.044E-03 + a2 * 3.687E-04 - a3 * 2.814E-06 + a4 * 8.938E-09;
      else
	attr[4] = 236.05828 - a * 2.81914E+00 + a2 * 8.39034E-03;
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);
      if (attr[0] > 179.0)
        sprintf(serr2, "magnitude value for Venus at phase angle i=%.1f is bad; formula is valid only for i < 179.0", attr[0]);
    } else if (ipl == SE_MARS) {
      double a = attr[0];
      double a2 = a * a;

      if (a <= 50.0)
	attr[4] = -1.601 + a * 0.02267 - a2 * 0.0001302;
      else
	attr[4] = -0.367 - a * 0.02573 + a2 * 0.0003445;
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);
    } else if (ipl == SE_JUPITER) {

      double a = attr[0];
      double a2 = a * a;
      attr[4] = -9.395 - a * 3.7E-04 + a2 * 6.16E-04;
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);
    } else if (ipl == SE_SATURN) {
      double a = attr[0];
      double sinB2;
      T = (tjd - dt - J2000) / 36525.0;
      in = (28.075216 - 0.012998 * T + 0.000004 * T * T) * DEGTORAD;
      om = (169.508470 + 1.394681 * T + 0.000412 * T * T) * DEGTORAD;

      sinB = (sin(in) * cos(lbr[1] * DEGTORAD)
                    * sin(lbr[0] * DEGTORAD - om)
                    - cos(in) * sin(lbr[1] * DEGTORAD));
      sinB2 = (sin(in) * cos(lbr2[1] * DEGTORAD)
                    * sin(lbr2[0] * DEGTORAD - om)
                    - cos(in) * sin(lbr2[1] * DEGTORAD));
      sinB = fabs(sin((asin(sinB) + asin(sinB2)) / 2.0));
      attr[4] = -8.914 - 1.825 * sinB + 0.026 * a - 0.378 * sinB * pow(2.7182818,-2.25 * a);
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);
    } else if (ipl == SE_URANUS) {

      double a = attr[0];
      double a2 = a * a;
      double fi_ = 0;
      attr[4] = -7.110 - 8.4E-04 * fi_ + a * 6.587E-3 + a2 * 1.045E-4;
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);

      attr[4] -= 0.05;
    } else if (ipl == SE_NEPTUNE) {
      if (tjd < 2444239.5) {
	attr[4] = -6.89;
      } else if (tjd <= 2451544.5) {
	attr[4] = -6.89 - 0.0055 * (tjd - 2444239.5) / 365.25;

      } else {
	attr[4] = -7.00;
      }
      attr[4] += 5 * log10(lbr2[2] * lbr[2]);
#else
    } else if (ipl == SE_SATURN) {
      double u1, u2, du;

      T = (tjd - dt - J2000) / 36525.0;
      in = (28.075216 - 0.012998 * T + 0.000004 * T * T) * DEGTORAD;
      om = (169.508470 + 1.394681 * T + 0.000412 * T * T) * DEGTORAD;
      sinB = fabs(sin(in) * cos(lbr[1] * DEGTORAD)
                    * sin(lbr[0] * DEGTORAD - om)
                    - cos(in) * sin(lbr[1] * DEGTORAD));
      u1 = atan2(sin(in) * tan(lbr2[1] * DEGTORAD)
                             + cos(in) * sin(lbr2[0] * DEGTORAD - om),
                        cos(lbr2[0] * DEGTORAD - om)) * RADTODEG;
      u2 = atan2(sin(in) * tan(lbr[1] * DEGTORAD)
                             + cos(in) * sin(lbr[0] * DEGTORAD - om),
                        cos(lbr[0] * DEGTORAD - om)) * RADTODEG;
      du = swe_degnorm(u1 - u2);
      if (du > 10)
        du = 360 - du;
      attr[4] = 5 * log10(lbr2[2] * lbr[2])
                  + mag_elem[ipl][1] * sinB
                  + mag_elem[ipl][2] * sinB * sinB
                  + mag_elem[ipl][3] * du
                  + mag_elem[ipl][0];
#endif
    } else if (ipl < SE_CHIRON) {
      attr[4] = 5 * log10(lbr2[2] * lbr[2])
                  + mag_elem[ipl][1] * attr[0] /100.0
                  + mag_elem[ipl][2] * attr[0] * attr[0] / 10000.0
                  + mag_elem[ipl][3] * attr[0] * attr[0] * attr[0] / 1000000.0
                  + mag_elem[ipl][0];
    } else if (ipl < NMAG_ELEM || ipl > SE_AST_OFFSET) {
      ph1 = pow(EULER, -3.33 * pow(tan(attr[0] * DEGTORAD / 2), 0.63));
      ph2 = pow(EULER, -1.87 * pow(tan(attr[0] * DEGTORAD / 2), 1.22));
      if (ipl < NMAG_ELEM) {
        me[0] = mag_elem[ipl][0];
        me[1] = mag_elem[ipl][1];
      } else if (ipl == SE_AST_OFFSET + 1566) {

                me[0] = 16.9;
                me[1] = 0.15;
      } else {
        me[0] = swed.ast_H;
        me[1] = swed.ast_G;
      }
      attr[4] = 5 * log10(lbr2[2] * lbr[2])
          + me[0]
          - 2.5 * log10((1 - me[1]) * ph1 + me[1] * ph2);
    } else {
      attr[4] = 0;
    }
  }
  if (ipl != SE_SUN && ipl != SE_EARTH) {

    if (swe_calc(tjd, SE_SUN, iflag | SEFLG_XYZ, xx2, serr) == ERR)
      return ERR;
    if (swe_calc(tjd, SE_SUN, iflag, lbr2, serr) == ERR)
      return ERR;
    attr[2] = acos(swi_dot_prod_unit(xx, xx2)) * RADTODEG;
  }

  if (ipl == SE_MOON) {
    double sinhp, xm[6];

    if (swe_calc(tjd, (int) ipl, epheflag|SEFLG_TRUEPOS|SEFLG_EQUATORIAL|SEFLG_RADIANS, xm, serr) == ERR)

      return ERR;
    sinhp = EARTH_RADIUS / xm[2] / AUNIT;
    attr[5] = asin(sinhp) / DEGTORAD;

    if (iflag & SEFLG_TOPOCTR) {
      if (swe_calc(tjd, (int) ipl, epheflag|SEFLG_XYZ|SEFLG_TOPOCTR, xm, serr) == ERR)
	return ERR;
      if (swe_calc(tjd, (int) ipl, epheflag|SEFLG_XYZ, xx, serr) == ERR)
	return ERR;
      attr[5] = acos(swi_dot_prod_unit(xm, xx)) / DEGTORAD;
#if 0
      {

      double tsid, h, e, f = EARTH_OBLATENESS;
      double cosz, sinz, phi;

      tsid = swe_sidtime(tjd - swe_deltat_ex(tjd, iflag, serr)) * 15 + swed.topd.geolon;

      h = swe_degnorm(tsid - xm[0] / DEGTORAD);

      e = sqrt(f * (2 - f));
      phi = atan((1 - e * e) * tan(swed.topd.geolat * DEGTORAD));

      cosz = sin(xm[1]) * sin(phi) + cos(xm[1]) * cos(phi) * cos(h * DEGTORAD);
      sinz = sqrt(1 - cosz * cosz);
      attr[5] = asin(sinz * sinhp / (1 - sinz * sinhp)) / DEGTORAD;
      }
#endif
    }
  }
  if (*serr2 != '\0' && serr != NULL)
    strcpy(serr, serr2);
  return iflag;
}

int32 CALL_CONV swe_pheno_ut(double tjd_ut, int32 ipl, int32 iflag, double *attr, char *serr)
{
  double deltat;
  int32 retflag = OK;
  int32 epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag == 0) {
    epheflag = SEFLG_SWIEPH;
    iflag |= SEFLG_SWIEPH;
  }
  deltat = swe_deltat_ex(tjd_ut, iflag, serr);
  retflag = swe_pheno(tjd_ut + deltat, ipl, iflag, attr, serr);

  if ((retflag & SEFLG_EPHMASK) != epheflag) {
    deltat = swe_deltat_ex(tjd_ut, retflag, serr);
    retflag = swe_pheno(tjd_ut + deltat, ipl, iflag, attr, serr);
  }
  return retflag;
}

static int find_maximum(double y00, double y11, double y2, double dx,
                        double *dxret, double *yret)
{
  double a, b, c, x, y;
  c = y11;
  b = (y2 - y00) / 2.0;
  a = (y2 + y00) / 2.0 - c;
  x = -b / 2 / a;
  y = (4 * a * c - b * b) / 4 / a;
  *dxret = (x - 1) * dx;
  if (yret != NULL)
    *yret = y;
  return OK;
}

static int find_zero(double y00, double y11, double y2, double dx,
                        double *dxret, double *dxret2)
{
  double a, b, c, x1, x2;
  c = y11;
  b = (y2 - y00) / 2.0;
  a = (y2 + y00) / 2.0 - c;
  if (b * b - 4 * a * c < 0)
    return ERR;
  x1 = (-b + sqrt(b * b - 4 * a * c)) / 2 / a;
  x2 = (-b - sqrt(b * b - 4 * a * c)) / 2 / a;
  *dxret = (x1 - 1) * dx;
  *dxret2 = (x2 - 1) * dx;
  return OK;
}

double rdi_twilight(int32 rsmi)
{
  double rdi = 0;
  if (rsmi & SE_BIT_CIVIL_TWILIGHT)
    rdi = 6;
  if (rsmi & SE_BIT_NAUTIC_TWILIGHT)
    rdi = 12;
  if (rsmi & SE_BIT_ASTRO_TWILIGHT)
    rdi = 18;
  return rdi;
}

static double get_sun_rad_plus_refr(int32 ipl, double dd, int32 rsmi, double refr)
{
  double rdi = 0;
  if (rsmi & SE_BIT_FIXED_DISC_SIZE) {
    if (ipl == SE_SUN)
      dd = 1.0;
    else if (ipl == SE_MOON)
      dd = 0.00257;
  }

  if (!(rsmi & SE_BIT_DISC_CENTER))
    rdi = asin( pla_diam[ipl] / 2.0 / AUNIT / dd) * RADTODEG;
  if (rsmi & SE_BIT_DISC_BOTTOM)
    rdi = -rdi;
  if (!(rsmi & SE_BIT_NO_REFRACTION)) {
    rdi += refr;
  }
  return rdi;
}

static int32 rise_set_fast(
               double tjd_ut, int32 ipl,
	       int32 epheflag, int32 rsmi,
               double *dgeo,
	       double atpress, double attemp,
               double *tret,
               char *serr)
{
  int i;
  double xx[6], xaz[6], xaz2[6];
  double dd, dt, refr;
  int32 iflag = epheflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  int32 iflagtopo = iflag | SEFLG_EQUATORIAL;
  double sda, armc, md, dmd, mdrise, rdi, tr, dalt;
  double decl;
  double tjd_ut0 = tjd_ut;
  int32 facrise = 1;
  int32 tohor_flag = SE_EQU2HOR;
  AS_BOOL is_second_run = FALSE;
  int nloop = 2;
  *tret = 0;
  if (ipl == SE_MOON)
    nloop = 4;
  if (rsmi & SE_CALC_SET)
    facrise = -1;
  if (!(rsmi & SE_BIT_GEOCTR_NO_ECL_LAT)) {
    iflagtopo |= SEFLG_TOPOCTR;
    swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);
  }
run_rise_again:
  if (swe_calc_ut(tjd_ut, ipl, iflagtopo, xx, serr) == ERR)
    return ERR;

  decl = xx[1];

  sda = -tan(dgeo[1] * DEGTORAD) * tan(decl * DEGTORAD);
  if (sda >= 1) {
    sda = 10;

  } else if (sda <= -1) {
    sda = 180;
  } else {
    sda = acos(sda) * RADTODEG;
  }

  armc = swe_degnorm(swe_sidtime(tjd_ut) * 15 + dgeo[0]);

  md = swe_degnorm(xx[0] - armc);
  mdrise = swe_degnorm(sda * facrise);

  dmd = swe_degnorm(md - mdrise);

#if 0
  if (dmd > 358) {
    tjd_ut -= 0.1;
    goto run_rise_again;
  }
#else
  if (dmd > 358) {
    dmd -= 360;
  }
#endif

  tr = tjd_ut + dmd / 360;

  rdi = 0;

  if (atpress == 0) {

    atpress = 1013.25 * pow(1 - 0.0065 * dgeo[2] / 288, 5.255);
  }
  swe_refrac_extended(0.000001, 0, atpress, attemp, const_lapse_rate, SE_APP_TO_TRUE, xx);
  refr = xx[1] - xx[0];

  if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT) {
    tohor_flag = SE_ECL2HOR;
    iflagtopo = iflag;
  } else {
    tohor_flag = SE_EQU2HOR;
    iflagtopo = iflag | SEFLG_EQUATORIAL;
    iflagtopo |= SEFLG_TOPOCTR;
    swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);
  }
  for (i = 0; i < nloop; i++) {
    if (swe_calc_ut(tr, ipl, iflagtopo, xx, serr) == ERR)
      return ERR;
    if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT)
      xx[1] = 0;
    rdi = get_sun_rad_plus_refr(ipl, xx[2], rsmi, refr);
    swe_azalt(tr, tohor_flag, dgeo, atpress, attemp, xx, xaz);
    swe_azalt(tr + 0.001, tohor_flag, dgeo, atpress, attemp, xx, xaz2);
    dd = (xaz2[1] - xaz[1]);
    dalt = xaz[1] + rdi;
    dt = dalt / dd / 1000.0;
    if (dt > 0.1) {
      dt = 0.1;
    } else if (dt < -0.1) {
      dt = -0.1;
    }
    if ((0) && fabs(dt) > 5.0 / 86400.0 && nloop < 20)
      nloop++;
    tr -= dt;
  }

  if (tr < tjd_ut0 && !is_second_run) {
    tjd_ut += 0.5;
    is_second_run = TRUE;
    goto run_rise_again;
  }
  *tret = tr;
  return OK;
}

#define SEFLG_EPHMASK	(SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)
int32 CALL_CONV swe_rise_trans(
               double tjd_ut, int32 ipl, char *starname,
	       int32 epheflag, int32 rsmi,
               double *geopos,
	       double atpress, double attemp,
               double *tret,
               char *serr)
{
  int32 retval = 0;

  AS_BOOL do_fixstar = (starname != NULL && *starname != '\0');
  if (!do_fixstar
    && (rsmi & (SE_CALC_RISE|SE_CALC_SET))
    && !(rsmi & SE_BIT_FORCE_SLOW_METHOD)
    && !(rsmi & (SE_BIT_CIVIL_TWILIGHT|SE_BIT_NAUTIC_TWILIGHT|SE_BIT_ASTRO_TWILIGHT))
    && (ipl >= SE_SUN && ipl <= SE_TRUE_NODE)
    && (fabs(geopos[1]) <= 60 || (ipl == SE_SUN && fabs(geopos[1]) <= 65))
    ) {
      retval = rise_set_fast(tjd_ut, ipl, epheflag, rsmi, geopos, atpress, attemp, tret, serr);
      return retval;
  }
  return swe_rise_trans_true_hor(tjd_ut, ipl, starname, epheflag, rsmi, geopos, atpress, attemp, 0, tret, serr);
}

int32 CALL_CONV swe_rise_trans_true_hor(
               double tjd_ut, int32 ipl, char *starname,
	       int32 epheflag, int32 rsmi,
               double *geopos,
	       double atpress, double attemp,
	       double horhgt,
               double *tret,
               char *serr)
{
  int i, j, k, ii, calc_culm, nculm = -1;
  double tjd_et = tjd_ut + swe_deltat_ex(tjd_ut, epheflag, serr);
  double xc[6], xh[20][6], ah[6], aha;
  double tculm[4], tcu, tc[20], h[20], t2[6], dc[6], dtint, dx, rdi, dd = 0;
  int32 iflag = epheflag;
  int jmax = 14;
  double t, te, tt, dt, twohrs = 1.0 / 12.0;
  double curdist;
  int32 tohor_flag = SE_EQU2HOR;
  int nazalt = 0;
  int ncalc = 0;
  AS_BOOL do_fixstar = (starname != NULL && *starname != '\0');
  if (geopos[2] < SEI_ECL_GEOALT_MIN || geopos[2] > SEI_ECL_GEOALT_MAX) {
    if (serr != NULL)
      sprintf(serr, "location for swe_rise_trans() must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }

  if (horhgt == -100) {
    horhgt = 0.0001 + calc_dip(geopos[2], atpress, attemp, const_lapse_rate);
  }

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  xh[0][0] = 0;

  iflag &= (SEFLG_EPHMASK | SEFLG_NONUT | SEFLG_TRUEPOS);
  *tret = 0;
  if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT) {
    tohor_flag = SE_ECL2HOR;
  } else {
    tohor_flag = SE_EQU2HOR;
    iflag |= SEFLG_EQUATORIAL;
    iflag |= SEFLG_TOPOCTR;
    swe_set_topo(geopos[0], geopos[1], geopos[2]);
  }
  if (rsmi & (SE_CALC_MTRANSIT | SE_CALC_ITRANSIT))
    return calc_mer_trans(tjd_ut, ipl, epheflag, rsmi,
		geopos, starname, tret, serr);
  if (!(rsmi & (SE_CALC_RISE | SE_CALC_SET)))
    rsmi |= SE_CALC_RISE;

  if (ipl == SE_SUN && (rsmi & (SE_BIT_CIVIL_TWILIGHT|SE_BIT_NAUTIC_TWILIGHT|SE_BIT_ASTRO_TWILIGHT))) {
    rsmi |= (SE_BIT_NO_REFRACTION | SE_BIT_DISC_CENTER);
    horhgt = -rdi_twilight(rsmi);

  }

  if (do_fixstar) {
    if (swe_fixstar(starname, tjd_et, iflag, xc, serr) == ERR)
      return ERR;
  }
  for (ii = 0, t = tjd_ut - twohrs; ii <= jmax; ii++, t += twohrs) {
    tc[ii] = t;
    if (!do_fixstar) {
      te = t + swe_deltat_ex(t, epheflag, serr);
      if (swe_calc(te, ipl, iflag, xc, serr) == ERR)
        return ERR;
ncalc++;
    }
    if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT)
      xc[1] = 0;

    if (ii == 0) {
      if (do_fixstar)
        dd = 0;
      else if (rsmi & SE_BIT_DISC_CENTER)
        dd = 0;
      else if (ipl < NDIAM)
        dd = pla_diam[ipl];
      else if (ipl > SE_AST_OFFSET)
        dd = swed.ast_diam * 1000;
      else
        dd = 0;
    }
    curdist = xc[2];
    if (rsmi & SE_BIT_FIXED_DISC_SIZE) {
      if (ipl == SE_SUN) {
        curdist = 1.0;
      } else if (ipl == SE_MOON) {
        curdist = 0.00257;
      }
    }

    rdi = asin( dd / 2 / AUNIT / curdist ) * RADTODEG;

    swe_azalt(t, tohor_flag, geopos, atpress, attemp, xc, xh[ii]);
nazalt++;
    if (rsmi & SE_BIT_DISC_BOTTOM) {

      xh[ii][1] -= rdi;
    } else {

      xh[ii][1] += rdi;
    }

    if (rsmi & SE_BIT_NO_REFRACTION) {
      xh[ii][1] -= horhgt;
      h[ii] = xh[ii][1];
    } else {
      swe_azalt_rev(t, SE_HOR2EQU, geopos, xh[ii], xc);
nazalt++;
      swe_azalt(t, SE_EQU2HOR, geopos, atpress, attemp, xc, xh[ii]);
nazalt++;
      xh[ii][1] -= horhgt;
      xh[ii][2] -= horhgt;
      h[ii] = xh[ii][2];
    }
    calc_culm = 0;
    if (ii > 1) {
      dc[0] = xh[ii-2][1];
      dc[1] = xh[ii-1][1];
      dc[2] = xh[ii][1];
      if (dc[1] > dc[0] && dc[1] > dc[2])
        calc_culm = 1;
      if (dc[1] < dc[0] && dc[1] < dc[2])
        calc_culm = 2;
    }
    if (calc_culm) {
      dt = twohrs;
      tcu = t - dt;
      find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dx);
      tcu += dtint + dt;
      dt /= 3;
      for (; dt > 0.0001; dt /= 3) {
        for (i = 0, tt = tcu - dt; i < 3; tt += dt, i++) {
          te = tt + swe_deltat_ex(tt, epheflag, serr);
          if (!do_fixstar) {
            if (swe_calc(te, ipl, iflag, xc, serr) == ERR)
              return ERR;
	  }
	  if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT)
	    xc[1] = 0;
	  ncalc++;
          swe_azalt(tt, tohor_flag, geopos, atpress, attemp, xc, ah);
	  nazalt++;
	  ah[1] -= horhgt;
          dc[i] = ah[1];
        }
        find_maximum(dc[0], dc[1], dc[2], dt, &dtint, &dx);
        tcu += dtint + dt;
      }
      nculm++;
      tculm[nculm] = tcu;
    }
  }

  for (i = 0; i <= nculm; i++) {
    for (j = 1; j <= jmax; j++) {
      if (tculm[i] < tc[j]) {
        for (k = jmax; k >= j; k--) {
          tc[k+1] = tc[k];
          h[k+1] = h[k];
        }
        tc[j] = tculm[i];
        if (!do_fixstar) {
          te = tc[j] + swe_deltat_ex(tc[j], epheflag, serr);
          if (swe_calc(te, ipl, iflag, xc, serr) == ERR)
            return ERR;
	  if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT)
	    xc[1] = 0;
ncalc++;
        }
        curdist = xc[2];
        if (rsmi & SE_BIT_FIXED_DISC_SIZE) {
          if ( ipl == SE_SUN ) {
            curdist = 1.0;
          } else if (ipl == SE_MOON) {
            curdist = 0.00257;
          }
        }

        rdi = asin( dd / 2 / AUNIT / curdist ) * RADTODEG;

        swe_azalt(tc[j], tohor_flag, geopos, atpress, attemp, xc, ah);
nazalt++;
        if (rsmi & SE_BIT_DISC_BOTTOM) {

          ah[1] -= rdi;
        } else {

	  ah[1] += rdi;
        }

	if (rsmi & SE_BIT_NO_REFRACTION) {
	  ah[1] -= horhgt;
	  h[j] = ah[1];
	} else {
	  swe_azalt_rev(tc[j], SE_HOR2EQU, geopos, ah, xc);
nazalt++;
	  swe_azalt(tc[j], SE_EQU2HOR, geopos, atpress, attemp, xc, ah);
nazalt++;
	  ah[1] -= horhgt;
	  ah[2] -= horhgt;
	  h[j] = ah[2];
	}
        jmax++;
        break;
      }
    }
  }
  *tret = 0;

  for (ii = 1; ii <= jmax; ii++) {
    if (h[ii-1] * h[ii] >= 0)
      continue;
    if (h[ii-1] < h[ii] && !(rsmi & SE_CALC_RISE))
      continue;
    if (h[ii-1] > h[ii] && !(rsmi & SE_CALC_SET))
      continue;
    dc[0] = h[ii-1];
    dc[1] = h[ii];
    t2[0] = tc[ii-1];
    t2[1] = tc[ii];
    for (i = 0; i < 20; i++) {
      t = (t2[0] + t2[1]) / 2;
      if (!do_fixstar) {
        te = t + swe_deltat_ex(t, epheflag, serr);
        if (swe_calc(te, ipl, iflag, xc, serr) == ERR)
          return ERR;
	if (rsmi & SE_BIT_GEOCTR_NO_ECL_LAT)
	  xc[1] = 0;
ncalc++;
      }
      curdist = xc[2];
      if (rsmi & SE_BIT_FIXED_DISC_SIZE) {
        if (ipl == SE_SUN) {
          curdist = 1.0;
        } else if (ipl == SE_MOON) {
          curdist = 0.00257;
        }
      }

      rdi = asin( dd / 2 / AUNIT / curdist ) * RADTODEG;

      swe_azalt(t, tohor_flag, geopos, atpress, attemp, xc, ah);
nazalt++;
      if (rsmi & SE_BIT_DISC_BOTTOM) {

        ah[1] -= rdi;
      } else {

	ah[1] += rdi;
      }

      if (rsmi & SE_BIT_NO_REFRACTION) {
	ah[1] -= horhgt;
	aha = ah[1];
      } else {
	swe_azalt_rev(t, SE_HOR2EQU, geopos, ah, xc);
	nazalt++;
	swe_azalt(t, SE_EQU2HOR, geopos, atpress, attemp, xc, ah);
	nazalt++;
	ah[1] -= horhgt;
	ah[2] -= horhgt;
	aha = ah[2];
      }
      if (aha * dc[0] <= 0) {
        dc[1] = aha;
        t2[1] = t;
      } else {
        dc[0] = aha;
        t2[0] = t;
      }
    }
    if (t > tjd_ut) {
     *tret = t;

     return OK;
    }
  }
  if (serr)
    sprintf(serr, "rise or set not found for planet %d", ipl);
  return -2;
}

static int32 calc_mer_trans(
               double tjd_ut, int32 ipl, int32 epheflag, int32 rsmi,
               double *geopos,
               char *starname,
               double *tret,
               char *serr)
{
  int i;
  double tjd_et = tjd_ut + swe_deltat_ex(tjd_ut, epheflag, serr);
  double armc, armc0, arxc, x0[6], x[6], t, te;
  double mdd;
  int32 iflag = epheflag;
  AS_BOOL do_fixstar = (starname != NULL && *starname != '\0');
  iflag &= SEFLG_EPHMASK;
  *tret = 0;
  iflag |= (SEFLG_EQUATORIAL | SEFLG_TOPOCTR);
  armc0 = swe_sidtime(tjd_ut) + geopos[0] / 15;
  if (armc0 >= 24)
    armc0 -= 24;
  if (armc0 < 0)
    armc0 += 24;
  armc0 *= 15;
  if (do_fixstar) {
    if (swe_fixstar(starname, tjd_et, iflag, x0, serr) == ERR)
      return ERR;
  } else {
    if (swe_calc(tjd_et, ipl, iflag, x0, serr) == ERR)
      return ERR;
  }

  x[0] = x0[0];
  x[1] = x0[1];
  t = tjd_ut;
  arxc = armc0;
  if (rsmi & SE_CALC_ITRANSIT)
    arxc = swe_degnorm(arxc + 180);
  for (i = 0; i < 4; i++) {
    mdd = swe_degnorm(x[0] - arxc);
    if (i > 0 && mdd > 180)
      mdd -= 360;
    t += mdd / 361;
    armc = swe_sidtime(t) + geopos[0] / 15;
    if (armc >= 24)
      armc -= 24;
    if (armc < 0)
      armc += 24;
    armc *= 15;
      arxc = armc;
    if (rsmi & SE_CALC_ITRANSIT)
      arxc = swe_degnorm(arxc + 180);
    if (!do_fixstar) {
      te = t + swe_deltat_ex(t, epheflag, serr);
      if (swe_calc(te, ipl, iflag, x, serr) == ERR)
	return ERR;
    }
  }
  *tret = t;
  return OK;
}

static const double el_node[8][4] =
  {{ 48.330893,  1.1861890,  0.00017587,  0.000000211,},
  { 76.679920,  0.9011190,  0.00040665, -0.000000080,},
  {  0       ,  0        ,  0         ,  0          ,},
  { 49.558093,  0.7720923,  0.00001605,  0.000002325,},
  {100.464441,  1.0209550,  0.00040117,  0.000000569,},
  {113.665524,  0.8770970, -0.00012067, -0.000002380,},
  { 74.005947,  0.5211258,  0.00133982,  0.000018516,},
  {131.784057,  1.1022057,  0.00026006, -0.000000636,},
  };
static const double el_peri[8][4] =
  {{ 77.456119,  1.5564775,  0.00029589,  0.000000056,},
  {131.563707,  1.4022188, -0.00107337, -0.000005315,},
  {102.937348,  1.7195269,  0.00045962,  0.000000499,},
  {336.060234,  1.8410331,  0.00013515,  0.000000318,},
  { 14.331309,  1.6126668,  0.00103127, -0.000004569,},
  { 93.056787,  1.9637694,  0.00083757,  0.000004899,},
  {173.005159,  1.4863784,  0.00021450,  0.000000433,},
  { 48.123691,  1.4262677,  0.00037918, -0.000000003,},
  };
static const double el_incl[8][4] =
  {{  7.004986,  0.0018215, -0.00001809,  0.000000053,},
  {  3.394662,  0.0010037, -0.00000088, -0.000000007,},
  {  0,         0,          0,           0          ,},
  {  1.849726, -0.0006010,  0.00001276, -0.000000006,},
  {  1.303270, -0.0054966,  0.00000465, -0.000000004,},
  {  2.488878, -0.0037363, -0.00001516,  0.000000089,},
  {  0.773196,  0.0007744,  0.00003749, -0.000000092,},
  {  1.769952, -0.0093082, -0.00000708,  0.000000028,},
  };
static const double el_ecce[8][4] =
  {{  0.20563175,  0.000020406, -0.0000000284, -0.00000000017,},
  {  0.00677188, -0.000047766,  0.0000000975,  0.00000000044,},
  {  0.01670862, -0.000042037, -0.0000001236,  0.00000000004,},
  {  0.09340062,  0.000090483, -0.0000000806, -0.00000000035,},
  {  0.04849485,  0.000163244, -0.0000004719, -0.00000000197,},
  {  0.05550862, -0.000346818, -0.0000006456,  0.00000000338,},
  {  0.04629590, -0.000027337,  0.0000000790,  0.00000000025,},
  {  0.00898809,  0.000006408, -0.0000000008, -0.00000000005,},
  };
static const double el_sema[8][4] =
  {{  0.387098310,  0.0,  0.0,  0.0,},
  {  0.723329820,  0.0,  0.0,  0.0,},
  {  1.000001018,  0.0,  0.0,  0.0,},
  {  1.523679342,  0.0,  0.0,  0.0,},
  {  5.202603191,  0.0000001913,  0.0,  0.0,},
  {  9.554909596,  0.0000021389,  0.0,  0.0,},
  { 19.218446062, -0.0000000372,  0.00000000098,  0.0,},
  { 30.110386869, -0.0000001663,  0.00000000069,  0.0,},
  };

static const double plmass[9] = {
    6023600,
     408523.719,
     328900.5,
    3098703.59,
       1047.348644,
       3497.9018,
      22902.98,
      19412.26,
  136566000,
};
static const int ipl_to_elem[15] = {2, 0, 0, 1, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 2,};
int32 CALL_CONV swe_nod_aps(double tjd_et, int32 ipl, int32 iflag,
                      int32  method,
                      double *xnasc, double *xndsc,
                      double *xperi, double *xaphe,
                      char *serr)
{
  int ij, i, j;
  int32 iplx;
  int32 ipli;
  int istart, iend;
  int32 iflJ2000;
  double daya, plm;
  double t = (tjd_et - J2000) / 36525, dt;
  double x[6], xx[24], *xp, xobs[6], x2000[6];
  double xpos[3][6], xnorm[6];
  double xposm[6];
  double xn[3][6], xs[3][6];
  double xq[3][6], xa[3][6];
  double xobs2[6], x2[6];
  double *xna, *xnd, *xpe, *xap;
  double incl, sema, ecce, parg, ea, vincl, vsema, vecce, pargx, eax;
  struct plan_data *pedp = &swed.pldat[SEI_EARTH];
  struct plan_data *psbdp = &swed.pldat[SEI_SUNBARY];
  struct plan_data pldat;
  double *xsun = psbdp->x;
  double *xear = pedp->x;
  const double *ep;
  double Gmsm, dzmin;
  double rxy, rxyz, fac, sgn;
  double sinnode, cosnode, sinincl, cosincl, sinu, cosu, sinE, cosE, cosE2;
  double uu, ny, ny2, c2, v2, pp, ro, ro2, rn, rn2;
  struct epsilon *oe;
  AS_BOOL is_true_nodaps = FALSE;
  AS_BOOL do_aberr = !(iflag & (SEFLG_TRUEPOS | SEFLG_NOABERR));
  AS_BOOL do_defl = !(iflag & SEFLG_TRUEPOS) && !(iflag & SEFLG_NOGDEFL);
  AS_BOOL do_focal_point = method & SE_NODBIT_FOPOINT;
  AS_BOOL ellipse_is_bary = FALSE;
  int32 iflg0;
  iflag &= ~(SEFLG_JPLHOR | SEFLG_JPLHOR_APPROX);

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  xna = xx;
  xnd = xx+6;
  xpe = xx+12;
  xap = xx+18;
  xpos[0][0] = 0;

  swi_force_app_pos_etc();
  method %= SE_NODBIT_FOPOINT;
  ipli = ipl;
  if (ipl == SE_SUN)
    ipli = SE_EARTH;
  if (ipl == SE_MOON) {
    do_defl = FALSE;
    if (!(iflag & SEFLG_HELCTR))
      do_aberr = FALSE;
  }
  iflg0 = (iflag & (SEFLG_EPHMASK|SEFLG_NONUT)) | SEFLG_SPEED | SEFLG_TRUEPOS;
  if (ipli != SE_MOON)
    iflg0 |= SEFLG_HELCTR;
  if (ipl == SE_MEAN_NODE || ipl == SE_TRUE_NODE ||
	  ipl == SE_MEAN_APOG || ipl == SE_OSCU_APOG ||
	  ipl < 0 ||
	  (ipl >= SE_NPLANETS && ipl <= SE_AST_OFFSET))

	  {
    if (serr != NULL)
      sprintf(serr, "nodes/apsides for planet %5.0f are not implemented", (double) ipl);
    if (xnasc != NULL)
      for (i = 0; i <= 5; i++)
	xnasc[i] = 0;
    if (xndsc != NULL)
      for (i = 0; i <= 5; i++)
	xndsc[i] = 0;
    if (xaphe != NULL)
      for (i = 0; i <= 5; i++)
	xaphe[i] = 0;
    if (xperi != NULL)
      for (i = 0; i <= 5; i++)
	xperi[i] = 0;
    return ERR;
  }
  for (i = 0; i < 24; i++)
    xx[i] = 0;

  if ((method == 0 || (method & SE_NODBIT_MEAN)) &&
        ((ipl >= SE_SUN && ipl <= SE_NEPTUNE) || ipl == SE_EARTH)) {
    if (ipl == SE_MOON) {
      swi_mean_lunar_elements(tjd_et, &xna[0], &xna[3], &xpe[0], &xpe[3]);
      incl = MOON_MEAN_INCL;
      vincl = 0;
      ecce = MOON_MEAN_ECC;
      vecce = 0;
      sema = MOON_MEAN_DIST / AUNIT;
      vsema = 0;
    } else {
      iplx = ipl_to_elem[ipl];
      ep = el_incl[iplx];
      incl = ep[0] + ep[1] * t + ep[2] * t * t + ep[3] * t * t * t;
      vincl = ep[1] / 36525;
      ep = el_sema[iplx];
      sema = ep[0] + ep[1] * t + ep[2] * t * t + ep[3] * t * t * t;
      vsema = ep[1] / 36525;
      ep = el_ecce[iplx];
      ecce = ep[0] + ep[1] * t + ep[2] * t * t + ep[3] * t * t * t;
      vecce = ep[1] / 36525;
      ep = el_node[iplx];

      xna[0] = ep[0] + ep[1] * t + ep[2] * t * t + ep[3] * t * t * t;
      xna[3] = ep[1] / 36525;

      ep = el_peri[iplx];
      xpe[0] = ep[0] + ep[1] * t + ep[2] * t * t + ep[3] * t * t * t;
      xpe[3] = ep[1] / 36525;
    }

    xnd[0] = swe_degnorm(xna[0] + 180);
    xnd[3] = xna[3];

    parg = xpe[0] = swe_degnorm(xpe[0] - xna[0]);
    pargx = xpe[3] = swe_degnorm(xpe[0] + xpe[3]  - xna[3]);

    swe_cotrans(xpe, xpe, -incl);

    swe_cotrans(xpe+3, xpe+3, -incl-vincl);

    xpe[0] = swe_degnorm(xpe[0] + xna[0]);

    xpe[3] = swe_degnorm(xpe[3] + xna[0] + xna[3]);

    xpe[3] = swe_degnorm(xpe[3] - xpe[0]);

    xpe[2] = sema * (1 - ecce);
    xpe[5] = (sema + vsema) * (1 - ecce - vecce) - xpe[2];

    xap[0] = swe_degnorm(xpe[0] + 180);
    xap[1] = -xpe[1];
    xap[3] = xpe[3];
    xap[4] = -xpe[4];
    if (do_focal_point) {
      xap[2] = sema * ecce * 2;
      xap[5] = (sema + vsema) * (ecce + vecce) * 2 - xap[2];
    } else {
      xap[2] = sema * (1 + ecce);
      xap[5] = (sema + vsema) * (1 + ecce + vecce) - xap[2];
    }

    ea = atan(tan(-parg * DEGTORAD / 2) * sqrt((1-ecce)/(1+ecce))) * 2;
    eax = atan(tan(-pargx * DEGTORAD / 2) * sqrt((1-ecce-vecce)/(1+ecce+vecce))) * 2;
    xna[2] = sema * (cos(ea) - ecce) / cos(parg * DEGTORAD);
    xna[5] = (sema+vsema) * (cos(eax) - ecce - vecce) / cos(pargx * DEGTORAD);
    xna[5] -= xna[2];
    ea = atan(tan((180 - parg) * DEGTORAD / 2) * sqrt((1-ecce)/(1+ecce))) * 2;
    eax = atan(tan((180 - pargx) * DEGTORAD / 2) * sqrt((1-ecce-vecce)/(1+ecce+vecce))) * 2;
    xnd[2] = sema * (cos(ea) - ecce) / cos((180 - parg) * DEGTORAD);
    xnd[5] = (sema+vsema) * (cos(eax) - ecce - vecce) / cos((180 - pargx) * DEGTORAD);
    xnd[5] -= xnd[2];

    for (i = 0, xp = xx; i < 4; i++, xp += 6) {

      xp[0] *= DEGTORAD;
      xp[1] *= DEGTORAD;
      xp[3] *= DEGTORAD;
      xp[4] *= DEGTORAD;
      swi_polcart_sp(xp, xp);
    }

  } else {

    if (swe_calc(tjd_et, ipli, iflg0, x, serr) == ERR)
      return ERR;
    iflJ2000 = (iflag & SEFLG_EPHMASK)|SEFLG_J2000|SEFLG_EQUATORIAL|SEFLG_XYZ|SEFLG_TRUEPOS|SEFLG_NONUT|SEFLG_SPEED;
    ellipse_is_bary = FALSE;
    if (ipli != SE_MOON) {
      if ((method & SE_NODBIT_OSCU_BAR) && x[2] > 6) {
        iflJ2000 |= SEFLG_BARYCTR;
        ellipse_is_bary = TRUE;
      } else {
        iflJ2000 |= SEFLG_HELCTR;
      }
    }

    if (ipli == SE_MOON) {
      dt = NODE_CALC_INTV;
      dzmin = 1e-15;
      Gmsm = GEOGCONST * (1 + 1 / EARTH_MOON_MRAT) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;
    } else {
      if ((ipli >= SE_MERCURY && ipli <= SE_PLUTO) || ipli == SE_EARTH)
        plm = 1 / plmass[ipl_to_elem[ipl]];
      else
        plm = 0;
      dt = NODE_CALC_INTV * 10 * x[2];
      dzmin = 1e-15 * dt / NODE_CALC_INTV;
      Gmsm = HELGRAVCONST * (1 + plm) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;
    }
    if (iflag & SEFLG_SPEED) {
      istart = 0;
      iend = 2;
    } else {
      istart = iend = 0;
      dt = 0;
    }
    for (i = istart, t = tjd_et - dt; i <= iend; i++, t += dt) {
      if (istart == iend)
        t = tjd_et;
      if (swe_calc(t, ipli, iflJ2000, xpos[i], serr) == ERR)
        return ERR;

      if (ipli == SE_EARTH) {
        if (swe_calc(t, SE_MOON, iflJ2000 & ~(SEFLG_BARYCTR|SEFLG_HELCTR), xposm, serr) == ERR)
          return ERR;
        for (j = 0; j <= 5; j++)
          xpos[i][j] += xposm[j] / (EARTH_MOON_MRAT + 1.0);
      }
      swi_plan_for_osc_elem(iflg0, t, xpos[i]);
    }
    for (i = istart; i <= iend; i++) {
      if (fabs(xpos[i][5]) < dzmin)
        xpos[i][5] = dzmin;
      fac = xpos[i][2] / xpos[i][5];
      sgn = xpos[i][5] / fabs(xpos[i][5]);
      for (j = 0; j <= 2; j++) {
        xn[i][j] = (xpos[i][j] - fac * xpos[i][j+3]) * sgn;
        xs[i][j] = -xn[i][j];
      }
    }
    for (i = istart; i <= iend; i++) {

      rxy =  sqrt(xn[i][0] * xn[i][0] + xn[i][1] * xn[i][1]);
      cosnode = xn[i][0] / rxy;
      sinnode = xn[i][1] / rxy;

      swi_cross_prod(xpos[i], xpos[i]+3, xnorm);
      rxy =  xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
      c2 = (rxy + xnorm[2] * xnorm[2]);
      rxyz = sqrt(c2);
      rxy = sqrt(rxy);
      sinincl = rxy / rxyz;
      cosincl = sqrt(1 - sinincl * sinincl);
      if (xnorm[2] < 0) cosincl = -cosincl;

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

      xq[i][0] = swi_mod2PI(uu - ny);
      xq[i][1] = 0;
      xq[i][2] = sema * (1 - ecce);

      swi_polcart(xq[i], xq[i]);
      swi_coortrf2(xq[i], xq[i], -sinincl, cosincl);
      swi_cartpol(xq[i], xq[i]);

      xq[i][0] += atan2(sinnode, cosnode);
      xa[i][0] = swi_mod2PI(xq[i][0] + PI);
      xa[i][1] = -xq[i][1];
      if (do_focal_point) {
        xa[i][2] = sema * ecce * 2;
      } else {
        xa[i][2] = sema * (1 + ecce);
      }
      swi_polcart(xq[i], xq[i]);
      swi_polcart(xa[i], xa[i]);

      ny = swi_mod2PI(ny - uu);
      ny2 = swi_mod2PI(ny + PI);

      cosE = cos(2 * atan(tan(ny / 2) / sqrt((1+ecce) / (1-ecce))));
      cosE2 = cos(2 * atan(tan(ny2 / 2) / sqrt((1+ecce) / (1-ecce))));

      rn = sema * (1 - ecce * cosE);
      rn2 = sema * (1 - ecce * cosE2);

      ro = sqrt(square_sum(xn[i]));
      ro2 = sqrt(square_sum(xs[i]));

      for (j = 0; j <= 2; j++) {
        xn[i][j] *= rn / ro;
        xs[i][j] *= rn2 / ro2;
      }
    }
    for (i = 0; i <= 2; i++) {
      if (iflag & SEFLG_SPEED) {
        xpe[i] = xq[1][i];
        xpe[i+3] = (xq[2][i] - xq[0][i]) / dt / 2;
        xap[i] = xa[1][i];
        xap[i+3] = (xa[2][i] - xa[0][i]) / dt / 2;
        xna[i] = xn[1][i];
        xna[i+3] = (xn[2][i] - xn[0][i]) / dt / 2;
        xnd[i] = xs[1][i];
        xnd[i+3] = (xs[2][i] - xs[0][i]) / dt / 2;
      } else {
        xpe[i] = xq[0][i];
        xpe[i+3] = 0;
        xap[i] = xa[0][i];
        xap[i+3] = 0;
        xna[i] = xn[0][i];
        xna[i+3] = 0;
        xnd[i] = xs[0][i];
        xnd[i+3] = 0;
      }
    }
    is_true_nodaps = TRUE;
  }

  if (ipli == SE_MOON && (iflag & (SEFLG_HELCTR | SEFLG_BARYCTR))) {
    swi_force_app_pos_etc();
    if (swe_calc(tjd_et, SE_SUN, iflg0, x, serr) == ERR)
      return ERR;
  } else {
    if (swe_calc(tjd_et, ipli, iflg0 | (iflag & SEFLG_TOPOCTR), x, serr) == ERR)
      return ERR;
  }

  if (iflag & SEFLG_TOPOCTR) {

    if (swi_get_observer(tjd_et, iflag, FALSE, xobs, serr) != OK)
      return ERR;

  } else {
    for (i = 0; i <= 5; i++)
      xobs[i] = 0;
  }
  if (iflag & (SEFLG_HELCTR | SEFLG_BARYCTR)) {
    if ((iflag & SEFLG_HELCTR) && !(iflag & SEFLG_MOSEPH))
      for (i = 0; i <= 5; i++)
        xobs[i] = xsun[i];
  } else if (ipl == SE_SUN && !(iflag & SEFLG_MOSEPH)) {
    for (i = 0; i <= 5; i++)
      xobs[i] = xsun[i];
  } else {

    for (i = 0; i <= 5; i++)
      xobs[i] += xear[i];
  }

  if (iflag & SEFLG_J2000)
    oe = &swed.oec2000;
  else
    oe = &swed.oec;

  for (ij = 0, xp = xx; ij < 4; ij++, xp += 6) {

    if (ipli == SE_EARTH && ij <= 1) {
      for (i = 0; i <= 5; i++)
    	xp[i] = 0;
      continue;
    }

    if (is_true_nodaps && !(iflag & SEFLG_NONUT)) {
      swi_coortrf2(xp, xp, -swed.nut.snut, swed.nut.cnut);
      if (iflag & SEFLG_SPEED)
        swi_coortrf2(xp+3, xp+3, -swed.nut.snut, swed.nut.cnut);
    }
    swi_coortrf2(xp, xp, -oe->seps, oe->ceps);
    swi_coortrf2(xp+3, xp+3, -oe->seps, oe->ceps);
    if (is_true_nodaps) {

      if (!(iflag & SEFLG_NONUT))
	swi_nutate(xp, iflag, TRUE);
    }

    swi_precess(xp, tjd_et, iflag, J_TO_J2000);
    if (iflag & SEFLG_SPEED)
      swi_precess_speed(xp, tjd_et, iflag, J_TO_J2000);

    if (ipli == SE_MOON) {
      for (i = 0; i <= 5; i++)
        xp[i] += xear[i];
    } else {
      if (!(iflag & SEFLG_MOSEPH) && !ellipse_is_bary)
        for (j = 0; j <= 5; j++)
          xp[j] += xsun[j];
    }

    for (j = 0; j <= 5; j++)
      xp[j] -= xobs[j];

    if (ipl == SE_SUN && !(iflag & (SEFLG_HELCTR | SEFLG_BARYCTR)))
      for (j = 0; j <= 5; j++)
        xp[j] = -xp[j];

    dt = sqrt(square_sum(xp)) * AUNIT / CLIGHT / 86400.0;
    if (do_defl)
      swi_deflect_light(xp, dt, iflag);

    if (do_aberr) {
      swi_aberr_light(xp, xobs, iflag);

      if (iflag & SEFLG_SPEED) {

        if (swe_calc(tjd_et - dt, ipli, iflg0 | (iflag & SEFLG_TOPOCTR), x2, serr) == ERR)
          return ERR;
        if (iflag & SEFLG_TOPOCTR) {

          for (i = 0; i <= 5; i++)
            xobs2[i] = swed.topd.xobs[i];
        } else {
          for (i = 0; i <= 5; i++)
            xobs2[i] = 0;
        }
        if (iflag & (SEFLG_HELCTR | SEFLG_BARYCTR)) {
          if ((iflag & SEFLG_HELCTR) && !(iflag & SEFLG_MOSEPH))
            for (i = 0; i <= 5; i++)
              xobs2[i] = xsun[i];
        } else if (ipl == SE_SUN && !(iflag & SEFLG_MOSEPH)) {
          for (i = 0; i <= 5; i++)
            xobs2[i] = xsun[i];
        } else {

          for (i = 0; i <= 5; i++)
            xobs2[i] += xear[i];
        }
        for (i = 3; i <= 5; i++)
          xp[i] += xobs[i] - xobs2[i];

        if (swe_calc(tjd_et, SE_SUN, iflg0 | (iflag & SEFLG_TOPOCTR), x2, serr) == ERR)
          return ERR;
      }
    }

    for (j = 0; j <= 5; j++)
      x2000[j] = xp[j];
    if (!(iflag & SEFLG_J2000)) {
      swi_precess(xp, tjd_et, iflag, J2000_TO_J);
      if (iflag & SEFLG_SPEED)
        swi_precess_speed(xp, tjd_et, iflag, J2000_TO_J);
    }

    if (!(iflag & SEFLG_NONUT))
      swi_nutate(xp, iflag, FALSE);

    for (j = 0; j <= 5; j++)
      pldat.xreturn[18+j] = xp[j];

    swi_coortrf2(xp, xp, oe->seps, oe->ceps);
    if (iflag & SEFLG_SPEED)
      swi_coortrf2(xp+3, xp+3, oe->seps, oe->ceps);
    if (!(iflag & SEFLG_NONUT)) {
      swi_coortrf2(xp, xp, swed.nut.snut, swed.nut.cnut);
      if (iflag & SEFLG_SPEED)
        swi_coortrf2(xp+3, xp+3, swed.nut.snut, swed.nut.cnut);
    }

    for (j = 0; j <= 5; j++)
      pldat.xreturn[6+j] = xp[j];

    if (iflag & SEFLG_SIDEREAL) {

      if (swed.sidd.sid_mode & SE_SIDBIT_ECL_T0) {
        if (swi_trop_ra2sid_lon(x2000, pldat.xreturn+6, pldat.xreturn+18, iflag) != OK)
          return ERR;

      } else if (swed.sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
        if (swi_trop_ra2sid_lon_sosy(x2000, pldat.xreturn+6, iflag) != OK)
          return ERR;
      } else {

        swi_cartpol_sp(pldat.xreturn+6, pldat.xreturn);
	if (swi_get_ayanamsa_ex(tjd_et, iflag, &daya, serr) == ERR)
	  return ERR;
        pldat.xreturn[0] -= daya * DEGTORAD;
        swi_polcart_sp(pldat.xreturn, pldat.xreturn+6);
      }
    }
    if ((iflag & SEFLG_XYZ) && (iflag & SEFLG_EQUATORIAL)) {
      for (j = 0; j <= 5; j++)
        xp[j] = pldat.xreturn[18+j];
      continue;
    }
    if (iflag & SEFLG_XYZ) {
      for (j = 0; j <= 5; j++)
        xp[j] = pldat.xreturn[6+j];
      continue;
    }

    swi_cartpol_sp(pldat.xreturn+18, pldat.xreturn+12);
    swi_cartpol_sp(pldat.xreturn+6, pldat.xreturn);

    if (!(iflag & SEFLG_RADIANS)) {
      for (j = 0; j < 2; j++) {
	pldat.xreturn[j] *= RADTODEG;
	pldat.xreturn[j+3] *= RADTODEG;
	pldat.xreturn[j+12] *= RADTODEG;
	pldat.xreturn[j+15] *= RADTODEG;
      }
    }
    if (iflag & SEFLG_EQUATORIAL) {
      for (j = 0; j <= 5; j++)
        xp[j] = pldat.xreturn[12+j];
      continue;
    } else {
      for (j = 0; j <= 5; j++)
        xp[j] = pldat.xreturn[j];
      continue;
    }
  }
  for (i = 0; i <= 5; i++) {
    if (i > 2 && !(iflag & SEFLG_SPEED))
      xna[i] = xnd[i] = xpe[i] = xap[i] = 0;
    if (xnasc != NULL)
      xnasc[i] = xna[i];
    if (xndsc != NULL)
      xndsc[i] = xnd[i];
    if (xperi != NULL)
      xperi[i] = xpe[i];
    if (xaphe != NULL)
      xaphe[i] = xap[i];
  }
  return OK;
}

int32 CALL_CONV swe_nod_aps_ut(double tjd_ut, int32 ipl, int32 iflag,
                      int32  method,
                      double *xnasc, double *xndsc,
                      double *xperi, double *xaphe,
                      char *serr) {

  return swe_nod_aps(tjd_ut + swe_deltat_ex(tjd_ut, iflag, serr),
                      ipl, iflag, method, xnasc, xndsc, xperi, xaphe,
                      serr);
}

#ifdef TEST_ORBEL_AA

static const double Gmsm_factor_AA[] = {
1,
0.9999999,
0.9999941,
0.9999941,
0.9999941,
0.9990404,
0.9987549,
0.998711,
0.99866025,
};
#endif
static int32 get_gmsm(double tjd_et, int32 ipl, int32 iflag, double r, double *gmsm, char *serr)
{
  int j;
  double Gmsm = 0, plm = 0, x[6];
  int32 iflJ2000p = (iflag & (SEFLG_EPHMASK |SEFLG_HELCTR|SEFLG_BARYCTR))|SEFLG_J2000|SEFLG_TRUEPOS|SEFLG_NONUT;
  if (!(iflJ2000p & (SEFLG_HELCTR|SEFLG_BARYCTR)))
    iflJ2000p |= SEFLG_HELCTR;
  if (ipl == SE_MOON) {
    Gmsm = GEOGCONST * (1 + 1 / EARTH_MOON_MRAT) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;
  } else {
    if ((ipl >= SE_MERCURY && ipl <= SE_PLUTO) || ipl == SE_EARTH) {
      plm = 0;

      if (iflag & SEFLG_ORBEL_AA) {
	if (ipl == SE_EARTH) {
	  plm = 1.0 / plmass[ipl_to_elem[ipl]];
	  plm += 1.0 / plmass[ipl_to_elem[SE_VENUS]];
	  plm += 1.0 / plmass[ipl_to_elem[SE_MERCURY]];
	} else {
	  for (j = ipl; j >= SE_MERCURY; j--) {
	    plm += 1.0 / plmass[ipl_to_elem[j]];
	  }
	  if (ipl >= SE_MARS)
	    plm += 1.0 / plmass[ipl_to_elem[SE_EARTH]];
	}

      } else {
	plm = 1.0 / plmass[ipl_to_elem[ipl]];
      }
      Gmsm = HELGRAVCONST * (1 + plm) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;
#ifdef TEST_ORBEL_AA
      if (!(iflag & SEFLG_ORBEL_AA))
	Gmsm /= Gmsm_factor_AA[ipl_to_elem[ipl]];
#endif

    } else {
      plm = 0;
      if (iflag & SEFLG_ORBEL_AA) {
	for (j = SE_MERCURY; j <= SE_PLUTO; j++) {
	  if (swe_calc(tjd_et, j, iflJ2000p, x, serr) == ERR)
	    return ERR;
	  if (r > x[2])
	    plm += 1.0 / plmass[ipl_to_elem[j]];
	}
	if (swe_calc(tjd_et, SE_EARTH, iflJ2000p, x, serr) == ERR)
	  return ERR;
	if (r > x[2])
	  plm += 1.0 / plmass[ipl_to_elem[SE_EARTH]];
      }
      Gmsm = HELGRAVCONST * (1 + plm) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;
    }
  }
  *gmsm = Gmsm;
  return OK;
}

int32 CALL_CONV swe_get_orbital_elements(
  double tjd_et,
  int32 ipl, int32 iflag,
  double *dret,
  char *serr)
{
  int j;
  double x[6], xpos[6], xposm[6], xn[6], xs[6], xnorm[6], xq[6], xa[6] ;

  int32 iflJ2000 = (iflag & SEFLG_EPHMASK)|SEFLG_J2000|SEFLG_XYZ|SEFLG_TRUEPOS|SEFLG_NONUT|SEFLG_SPEED;
  int32 iflJ2000p = (iflag & SEFLG_EPHMASK)|SEFLG_J2000|SEFLG_TRUEPOS|SEFLG_NONUT|SEFLG_SPEED;
  double Gmsm;

  double fac, sgn, rxy, rxyz, c2, cosnode, sinnode;
  double incl, node, parg, peri, mlon;
  double csid, ctro, csyn, dmot, pa;
  double ytrop, ysid, T, T2, T3, T4, T5;
  double sinincl, cosincl, cosu, sinu, uu, eanom, tanom, manom;
  double v2, sema, pp, ecce, cosE, sinE, ny, ny2, rn, rn2, ro, ro2, cosE2;
  double r, ecce2;
  if (ipl <= 0 || ipl == SE_MEAN_NODE || ipl == SE_TRUE_NODE || ipl == SE_MEAN_APOG || ipl == SE_OSCU_APOG || ipl == SE_INTP_APOG || ipl == SE_INTP_PERG) {
    if (serr != NULL)
      sprintf(serr, "error in swe_get_orbital_elements(): object %d not valid\n", ipl);
    return ERR;
  }

  if (swe_calc(tjd_et, ipl, iflJ2000p, x, serr) == ERR)
    return ERR;
  r =  x[2];
  if (ipl != SE_MOON) {
    if ((iflag & SEFLG_BARYCTR) && r > 6) {
      iflJ2000 |= SEFLG_BARYCTR;
    } else {
      iflJ2000 |= SEFLG_HELCTR;
    }
  }
  if (get_gmsm(tjd_et, ipl, iflag, r, &Gmsm, serr))
    return ERR;
  if (swe_calc(tjd_et, ipl, iflJ2000, xpos, serr) == ERR)
    return ERR;

  if (ipl == SE_EARTH) {
    if (swe_calc(tjd_et, SE_MOON, iflJ2000 & ~(SEFLG_BARYCTR|SEFLG_HELCTR), xposm, serr) == ERR)
      return ERR;
    for (j = 0; j <= 5; j++)
      xpos[j] += xposm[j] / (EARTH_MOON_MRAT + 1.0);
  }
  fac = xpos[2] / xpos[5];
  sgn = xpos[5] / fabs(xpos[5]);
  for (j = 0; j <= 2; j++) {
    xn[j] = (xpos[j] - fac * xpos[j+3]) * sgn;
    xs[j] = -xn[j];
  }

  rxy =  sqrt(xn[0] * xn[0] + xn[1] * xn[1]);
  cosnode = xn[0] / rxy;
  sinnode = xn[1] / rxy;

  swi_cross_prod(xpos, xpos+3, xnorm);
  rxy =  xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
  c2 = (rxy + xnorm[2] * xnorm[2]);
  rxyz = sqrt(c2);
  rxy = sqrt(rxy);
  sinincl = rxy / rxyz;
  cosincl = sqrt(1 - sinincl * sinincl);
  if (xnorm[2] < 0) cosincl = -cosincl;
  incl = acos(cosincl) * RADTODEG;

  cosu = xpos[0] * cosnode + xpos[1] * sinnode;
  sinu = xpos[2] / sinincl;
  uu = atan2(sinu, cosu);

  rxyz = sqrt(square_sum(xpos));
  v2 = square_sum((xpos+3));
  sema = 1.0 / (2.0 / rxyz - v2 / Gmsm);

  pp = c2 / Gmsm;
  ecce = pp / sema;
  if (ecce > 1)
    ecce = 1;
  ecce = sqrt(1 - ecce);

  ecce2 = ecce;
  if (ecce2 == 0)
    ecce2 = 0.0000000001;
  cosE = 1 / ecce2 * (1 - rxyz / sema);
  sinE = 1 / ecce2 / sqrt(sema * Gmsm) * dot_prod(xpos, (xpos+3));
  eanom = swe_degnorm(atan2(sinE, cosE) * RADTODEG);

  ny = 2 * atan(sqrt((1+ecce)/(1-ecce)) * sinE / (1 + cosE));
  tanom = swe_degnorm(ny * RADTODEG);
  if (eanom > 180 && tanom < 180)
    tanom += 180;
  if (eanom < 180 && tanom > 180)
    tanom -= 180;

  manom = swe_degnorm(eanom - ecce * RADTODEG * sin(eanom * DEGTORAD));

  xq[0] = swi_mod2PI(uu - ny);
  parg = xq[0] * RADTODEG;
  xq[1] = 0;
  xq[2] = sema * (1 - ecce);

  swi_polcart(xq, xq);
  swi_coortrf2(xq, xq, -sinincl, cosincl);
  swi_cartpol(xq, xq);

  xq[0] += atan2(sinnode, cosnode);
  xa[0] = swi_mod2PI(xq[0] + PI);
  xa[1] = -xq[1];

    xa[2] = sema * (1 + ecce);

  swi_polcart(xq, xq);
  swi_polcart(xa, xa);

  ny = swi_mod2PI(ny - uu);
  ny2 = swi_mod2PI(ny + PI);

  cosE = cos(2 * atan(tan(ny / 2) / sqrt((1+ecce) / (1-ecce))));
  cosE2 = cos(2 * atan(tan(ny2 / 2) / sqrt((1+ecce) / (1-ecce))));

  rn = sema * (1 - ecce * cosE);
  rn2 = sema * (1 - ecce * cosE2);

  ro = sqrt(square_sum(xn));
  ro2 = sqrt(square_sum(xs));

  for (j = 0; j <= 2; j++) {
    xn[j] *= rn / ro;
    xs[j] *= rn2 / ro2;
  }
  swi_cartpol(xn, xn);
  swi_cartpol(xq, xq);
  node = xn[0] * RADTODEG;
  peri = swe_degnorm(node + parg);
  mlon = swe_degnorm(manom + peri);
  csid = sema * sqrt(sema);
  if (ipl == SE_MOON) {
    double semam = sema * AUNIT / 383397772.5;
    csid = semam * sqrt(semam);
    csid *= 27.32166 / 365.25636300;
  }
  dmot = 0.9856076686 / csid;
  csid *= 365.25636 / 365.242189;

  T = (tjd_et - J2000) / 365250.0;
  T2 = T * T; T3 = T2 * T; T4 = T3 * T; T5 = T4 * T;
  pa = (50288.200 + 222.4045 * T + 0.2095 * T2 - 0.9408 * T3 - 0.0090 * T4 + 0.0010 * T5) / 3600.0 / 365250.0;

  ysid = (1295977422.83429 - 2 * 2.0441 * T  - 3 * 0.00523 * T * T) / 3600.0 / 365250.0;
  ysid = 360.0 / ysid;
  ytrop = (1296027711.03429 + 2 * 109.15809 * T  + 3 * 0.07207 * T2 - 4 * 0.23530 * T3 - 5 * 0.00180 * T4 + 6 * 0.00020 * T5) / 3600.0 / 365250.0;
  ytrop = 360.0 / ytrop;
  ctro = 360.0 / (dmot + pa) / 365.242189;
  ctro *= ysid / ytrop;
  if (ipl == SE_EARTH)
    csyn = 0;
  else
    csyn = 360.0 / (0.9856076686 - dmot);
  dret[0] = sema;
  dret[1] = ecce;
  dret[2] = incl;
  dret[3] = node;
  dret[4] = parg;
  dret[5] = peri;
  dret[6] = manom;
  dret[7] = tanom;
  dret[8] = eanom;
  dret[9] = mlon;
  dret[10] = csid;
  dret[11] = dmot;
  dret[12] = ctro;
  dret[13] = csyn;
  dret[14] = tjd_et - dret[6] / dmot;
  dret[15] = sema * (1 - ecce);
  dret[16] = sema * (1 + ecce);

  return OK;
}

static void osc_get_orbit_constants(double *dp, double *pqr)
{
  double sema = dp[0];
  double ecce = dp[1];
  double incl = dp[2];
  double node = dp[3];
  double parg = dp[4];
  double cosnode = cos(node * DEGTORAD);
  double sinnode = sin(node * DEGTORAD);
  double cosincl = cos(incl * DEGTORAD);
  double sinincl = sin(incl * DEGTORAD);
  double cosparg = cos(parg * DEGTORAD);
  double sinparg = sin(parg * DEGTORAD);
  double fac = sqrt((1 - ecce) * (1 + ecce));
  pqr[0] = cosparg * cosnode - sinparg * cosincl * sinnode;
  pqr[1] = -sinparg * cosnode - cosparg * cosincl * sinnode;
  pqr[2] = sinincl * sinnode;
  pqr[3] = cosparg * sinnode + sinparg * cosincl * cosnode;
  pqr[4] = -sinparg * sinnode + cosparg * cosincl * cosnode;
  pqr[5] = -sinincl * cosnode;
  pqr[6] = sinparg * sinincl;
  pqr[7] = cosparg * sinincl;
  pqr[8] = cosincl;
  pqr[9] = sema;
  pqr[10] = ecce;
  pqr[11] = fac;
}

static void osc_get_ecl_pos(double ean, double *pqr, double *xp)
{
  double x[2];
  double cose = cos(ean * DEGTORAD);
  double sine = sin(ean * DEGTORAD);
  double sema = pqr[9];
  double ecce = pqr[10];
  double fac = pqr[11];
  x[0] = sema * (cose - ecce);
  x[1] = sema * fac * sine;

  xp[0] = pqr[0] * x[0] + pqr[1] * x[1];
  xp[1] = pqr[3] * x[0] + pqr[4] * x[1];
  xp[2] = pqr[6] * x[0] + pqr[7] * x[1];
}

static double get_dist_from_2_vectors(double *x1, double *x2)
{
  double r0, r1, r2;
  r0 = x1[0] - x2[0];
  r1 = x1[1] - x2[1];
  r2 = x1[2] - x2[2];
  return sqrt(r0 * r0 + r1 * r1 + r2 * r2);
}

static void osc_iterate_max_dist(double ean, double *pqr, double *xa, double *xb, double *deanopt, double *drmax, AS_BOOL high_prec)
{
  int i;
  double r, rmax, eansv = 0, dstep, dstep_min = 1;
  if (high_prec)
    dstep_min = 0.000001;
  ean = 0;
  osc_get_ecl_pos(ean, pqr, xa);
  r = get_dist_from_2_vectors(xb, xa);
  rmax = r;
  dstep = 1;
  while (dstep >= dstep_min) {

    for (i = 0; i < 2; i++) {
      while(r >= rmax) {
	eansv = ean;
	if (i == 0)
	  ean += dstep;
	else
	  ean -= dstep;
	osc_get_ecl_pos(ean, pqr, xa);
	r = get_dist_from_2_vectors(xb, xa);
	if (r > rmax)
	  rmax = r;
      }
      ean = eansv;
      r = rmax;
    }
    ean = eansv;
    r = rmax;
    dstep /= 10;
  }
  *drmax = rmax;
  *deanopt = eansv;
}

static void osc_iterate_min_dist(double ean, double *pqr, double *xa, double *xb, double *deanopt, double *drmin, AS_BOOL high_prec)
{
  int i;
  double r, rmin, eansv = 0, dstep, dstep_min = 1;
  if (high_prec)
    dstep_min = 0.000001;
  ean = 0;
  osc_get_ecl_pos(ean, pqr, xa);
  r = get_dist_from_2_vectors(xb, xa);
  rmin = r;
  dstep = 1;
  while (dstep >= dstep_min) {

    for (i = 0; i < 2; i++) {
      while(r <= rmin) {
	eansv = ean;
	if (i == 0)
	  ean += dstep;
	else
	  ean -= dstep;
	osc_get_ecl_pos(ean, pqr, xa);
	r = get_dist_from_2_vectors(xb, xa);
	if (r < rmin)
	  rmin = r;
      }
      ean = eansv;
      r = rmin;
    }
    ean = eansv;
    r = rmin;
    dstep /= 10;
  }
  *drmin = rmin;
  *deanopt = eansv;
}

static int32 orbit_max_min_true_distance_helio(double tjd_et, int ipl, int32 iflag, double *dmax, double *dmin, double *dtrue, char *serr)
{
  double xinner[3], pqri[20];
  double eani;
  double de[50];
  int32 retval;
  int32 ipli = ipl;
  int32 iflagi = (iflag & (SEFLG_EPHMASK | SEFLG_HELCTR | SEFLG_BARYCTR));
  if (ipl == SE_SUN) {
    ipli = SE_EARTH;
  }

  if ((retval = swe_get_orbital_elements(tjd_et, ipli, iflagi, de, serr)) == ERR)
    return ERR;
  *dmax = de[16];
  *dmin = de[15];
  osc_get_orbit_constants(de, pqri);

  eani = de[8];

  osc_get_ecl_pos(eani, pqri, xinner);

  *dtrue = sqrt(xinner[0] * xinner[0] + xinner[1] * xinner[1] + xinner[2] * xinner[2]);
#ifdef DEBUG_REL_DIST
  printf("rtrue=%.17f (%.17f, %.17f\n", *dtrue, *dmin, *dmax);
#endif
  return retval;
}

int32 CALL_CONV swe_orbit_max_min_true_distance(double tjd_et, int32 ipl, int32 iflag, double *dmax, double *dmin, double *dtrue, char *serr)
{
  int i, j, k, retval;
  int32 iflagi = (iflag & (SEFLG_EPHMASK | SEFLG_HELCTR | SEFLG_BARYCTR));
  double dp[50], de[50];
  double xouter[3], xinner[3], max_xouter[3], max_xinner[3], min_xouter[3], min_xinner[3], pqro[20], pqri[20];
  double eano, eani;
  double *douter, *dinner;
  double r, rtrue, rmax = 0, rmin = 100000000, rminsv = 0, rmaxsv = 0;
  double min_eanisv = 0, min_eanosv = 0, max_eanisv = 0, max_eanosv = 0;
  int ncnt;
  double dstep;
  double nitermax = 300;

  if (ipl == SE_SUN || ipl == SE_MOON || (iflagi & (SEFLG_HELCTR | SEFLG_BARYCTR))) {
    retval = orbit_max_min_true_distance_helio(tjd_et, ipl, iflagi, dmax, dmin, dtrue, serr);
    return retval;
  }
  if ((retval = swe_get_orbital_elements(tjd_et, ipl, iflagi, dp, serr)) == ERR)
    return ERR;
  if ((retval = swe_get_orbital_elements(tjd_et, SE_EARTH, iflagi, de, serr)) == ERR)
    return ERR;
  if (de[0] > dp[0]) {
    douter = de;
    dinner = dp;
  } else {
    douter = dp;
    dinner = de;
  }
  osc_get_orbit_constants(douter, pqro);
  osc_get_orbit_constants(dinner, pqri);
  eano = douter[8];
  eani = dinner[8];
  osc_get_ecl_pos(eano, pqro, xouter);
  osc_get_ecl_pos(eani, pqri, xinner);
  rtrue = get_dist_from_2_vectors(xouter, xinner);

  ncnt = 182;
  dstep = 2;
  for (i = 0; i < 3; i++) {
    max_xouter[i] = 0;
    max_xinner[i] = 0;
    min_xouter[i] = 0;
    min_xinner[i] = 0;
  }
  for (j = 0; j < ncnt; j++) {
    eano = (double) j * dstep;
    osc_get_ecl_pos(eano, pqro, xouter);
    for (i = 0; i < ncnt; i++) {
      eani = (double) i;
      osc_get_ecl_pos(eani, pqri, xinner);
      r = get_dist_from_2_vectors(xouter, xinner);

      if (r > rmax) {
        rmax = r;
        max_eanisv = eani;
        max_eanosv = eano;
	for (k = 0; k < 3; k++) {
	  max_xouter[k] = xouter[k];
	  max_xinner[k] = xinner[k];
	}
      }
      if (r < rmin) {
        rmin = r;
        min_eanisv = eani;
        min_eanosv = eano;
	for (k = 0; k < 3; k++) {
	  min_xouter[k] = xouter[k];
	  min_xinner[k] = xinner[k];
	}
      }
    }
  }

  eani = max_eanisv;
  eano = max_eanosv;
  for (k = 0; k < 3; k++) {
    xouter[k] = max_xouter[k];
    xinner[k] = max_xinner[k];
  }
  for (k = 0; k <= nitermax; k++) {
    osc_iterate_max_dist(eani, pqri, xinner, xouter, &eani, &rmax, TRUE);
    osc_iterate_max_dist(eano, pqro, xouter, xinner, &eano, &rmax, TRUE);
    if (k > 0 && fabs(rmax - rmaxsv) < 0.00000001)
      break;
    rmaxsv = rmax;
  }

  eani = min_eanisv;
  eano = min_eanosv;
  for (k = 0; k < 3; k++) {
    xouter[k] = min_xouter[k];
    xinner[k] = min_xinner[k];
  }
  for (k = 0; k <= nitermax; k++) {
    osc_iterate_min_dist(eani, pqri, xinner, xouter, &eani, &rmin, TRUE);
    osc_iterate_min_dist(eano, pqro, xouter, xinner, &eano, &rmin, TRUE);
    if (k > 0 && fabs(rmin - rminsv) < 0.00000001)
      break;
    rminsv = rmin;
  }
  *dmax = rmax;
  *dmin = rmin;
  *dtrue = rtrue;
  return retval;
}

int32 CALL_CONV swe_gauquelin_sector(
  double t_ut,
  int32 ipl,

  char *starname,
  int32 iflag,
  int32 imeth,

  double *geopos,

  double atpress,

  double attemp,

  double *dgsect,
  char *serr)
{
  AS_BOOL rise_found = TRUE;
  AS_BOOL set_found = TRUE;
  int32 retval;
  double tret[3];
  double t_et, t;
  double x0[6];
  double eps, nutlo[2], armc;
  int32 epheflag = iflag & SEFLG_EPHMASK;
  AS_BOOL do_fixstar = (starname != NULL && *starname != '\0');
  int32 risemeth = 0;
  AS_BOOL above_horizon = FALSE;
  if (imeth < 0 || imeth > 5) {
    if (serr)
          sprintf(serr, "invalid method: %d", imeth);
    return ERR;
  }

  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;

  if (imeth == 0 || imeth == 1) {
    t_et = t_ut + swe_deltat_ex(t_ut, iflag, serr);
    eps = swi_epsiln(t_et, iflag) * RADTODEG;
    swi_nutation(t_et, iflag, nutlo);
    nutlo[0] *= RADTODEG;
    nutlo[1] *= RADTODEG;
    armc = swe_degnorm(swe_sidtime0(t_ut, eps + nutlo[1], nutlo[0]) * 15 + geopos[0]);
    if (do_fixstar) {
      if (swe_fixstar(starname, t_et, iflag, x0, serr) == ERR)
	return ERR;
    } else {
      if (swe_calc(t_et, ipl, iflag, x0, serr) == ERR)
	return ERR;
    }
    if (imeth == 1)
      x0[1] = 0;
    *dgsect = swe_house_pos(armc, geopos[1], eps + nutlo[1], 'G', x0, NULL);
    return OK;
  }

  if (imeth == 2 || imeth == 4)
    risemeth |= SE_BIT_NO_REFRACTION;
  if (imeth == 2 || imeth == 3)
    risemeth |= SE_BIT_DISC_CENTER;

  retval = swe_rise_trans(t_ut, ipl, starname, epheflag, SE_CALC_RISE|risemeth, geopos, atpress, attemp, &(tret[0]), serr);
  if (retval == ERR) {
    return ERR;
  } else if (retval == -2) {

    rise_found = FALSE;
  }

  retval = swe_rise_trans(t_ut, ipl, starname, epheflag, SE_CALC_SET|risemeth, geopos, atpress, attemp, &(tret[1]), serr);
  if (retval == ERR) {
    return ERR;
  } else if (retval == -2) {
    set_found = FALSE;
  }
  if (tret[0] < tret[1] && rise_found == TRUE) {
    above_horizon = FALSE;

    t = t_ut - 1.2;
    if (set_found) t = tret[1] - 1.2;
    set_found = TRUE;
    retval = swe_rise_trans(t, ipl, starname, epheflag, SE_CALC_SET|risemeth, geopos, atpress, attemp, &(tret[1]), serr);
    if (retval == ERR) {
      return ERR;
    } else if (retval == -2) {
      set_found = FALSE;
    }
  } else if (tret[0] >= tret[1] && set_found == TRUE) {
    above_horizon = TRUE;

    t = t_ut - 1.2;
    if (rise_found) t = tret[0] - 1.2;
    rise_found = TRUE;
    retval = swe_rise_trans(t, ipl, starname, epheflag, SE_CALC_RISE|risemeth, geopos, atpress, attemp, &(tret[0]), serr);
    if (retval == ERR) {
      return ERR;
    } else if (retval == -2) {
      rise_found = FALSE;
    }
  }
  if (rise_found && set_found) {
    if (above_horizon) {
      *dgsect = (t_ut - tret[0]) / (tret[1] - tret[0]) * 18 + 1;
    } else {
      *dgsect = (t_ut - tret[1]) / (tret[0] - tret[1]) * 18 + 19;
    }
    return OK;
  } else {
    *dgsect = 0;
    if (serr)
      sprintf(serr, "rise or set not found for planet %d", ipl);
    return ERR;
  }
}
