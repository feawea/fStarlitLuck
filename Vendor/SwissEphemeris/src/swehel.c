#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"
#include <sys/stat.h>

#define PLSV   0
#define criticalangle   0.0
#define BNIGHT   	1479.0
#define BNIGHT_FACTOR   1.0
#define PI		M_PI
#define Min2Deg   (1.0 / 60.0)
#define SWEHEL_DEBUG  0
#define DONE  1
#define MaxTryHours   4
#define TimeStepDefault	1
#define LocalMinStep	8

#define Y2D   365.25
#define D2Y   (1 / Y2D)
#define D2H   24.0
#define H2S   3600.0
#define D2S   (D2H * H2S)
#define S2H   (1.0 / H2S)
#define JC2D   36525.0
#define M2S   60.0

#define REFR_SINCLAIR    0
#define REFR_BENNETTH    1
#define FormAstroRefrac   REFR_SINCLAIR
#define GravitySource   2
#define REarthSource   1

#define StartYear   1820
#define Average   1.80546834626888
#define Periodicy   1443.67123144531
#define Amplitude   3.75606495492684
#define phase   0
#define MAX_COUNT_SYNPER           5
#define MAX_COUNT_SYNPER_MAX 1000000
#define AvgRadiusMoon  (15.541 / 60)

#define Ra   6378136.6
#define Rb   6356752.314

#define nL2erg 	(1.02E-15)
#define erg2nL 	(1 / nL2erg)
#define MoonDistance 	384410.4978
#define scaleHwater 	3000.0
#define scaleHrayleigh 	8515.0
#define scaleHaerosol 	3745.0
#define scaleHozone 	20000.0
#define astr2tau 	0.921034037197618
#define tau2astr 	1 / astr2tau

#define C2K   273.15
#define DELTA   18.36
#define TempNulDiff   0.000001
#define PressRef   1000
#define MD   28.964
#define MW   18.016
#define GCR   8314.472
#define LapseSA   0.0065
#define LapseDA   0.0098

#define LowestAppAlt   -3.5

#define epsilon   0.001

#define staticAirmass   0

#define GOpticMag   1
#define GOpticTrans   0.8
#define GBinocular   1
#define GOpticDia   50

static double mymin(double a, double b)
{
  if (a <= b)
    return a;
  return b;
}

static double mymax(double a, double b)
{
  if (a >= b)
    return a;
  return b;
}

static double Tanh(double x)
{
  return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
}

static double CVA(double B, double SN, int32 helflag)
{

  AS_BOOL is_scotopic = FALSE;

  if (B < 1394)
    is_scotopic = TRUE;
  if (helflag & SE_HELFLAG_VISLIM_PHOTOPIC)
    is_scotopic = FALSE;
  if (helflag & SE_HELFLAG_VISLIM_SCOTOPIC)
    is_scotopic = TRUE;
  if (is_scotopic)
    return mymin(900, 380 / SN * pow(10, (0.3 * pow(B, (-0.29))))) / 60.0 / 60.0;
  else
    return (40.0 / SN) * pow(10, (8.28 * pow(B, (-0.29)))) / 60.0 / 60.0;
}

static double PupilDia(double Age, double B)
{

  return (0.534 - 0.00211 * Age - (0.236 - 0.00127 * Age) * Tanh(0.4 * log(B) / log(10) - 2.2)) * 10;
}

static double OpticFactor(double Bback, double kX, double *dobs, double JDNDaysUT, char *ObjectName, int TypeFactor, int helflag)
{
  double Pst, CIb, CIi, ObjectSize, Fb, Fe, Fsc, Fci, Fcb, Ft, Fp, Fa, Fr, Fm;
  double Age = dobs[0];
  double SN = dobs[1], SNi;
  double Binocular = dobs[2];
  double OpticMag = dobs[3];
  double OpticDia = dobs[4];
  double OpticTrans = dobs[5];
  AS_BOOL is_scotopic = FALSE;
  JDNDaysUT += 0.0;
  SNi = SN;
  if (SNi <= 0.00000001) SNi = 0.00000001;

  Pst = PupilDia(23, Bback);
  if (OpticMag == 1) {
     OpticTrans = 1;
     OpticDia = Pst;
  }
#if 0
  if (OpticMag == 0) {
    OpticTrans = 1;
    OpticDia = Pst;
    Binocular = 1;
    OpticMag = 1;
  }
#endif

  CIb = 0.7;
  CIi = 0.5;
  ObjectSize = 0;
  if (strcmp(ObjectName, "moon") == 0) {

    ;
  }
  Fb = 1;
  if (Binocular == 0) Fb = 1.41;

  if (Bback < 1645)
    is_scotopic = TRUE;
  if (helflag & SE_HELFLAG_VISLIM_PHOTOPIC)
    is_scotopic = FALSE;
  if (helflag & SE_HELFLAG_VISLIM_SCOTOPIC)
    is_scotopic = TRUE;
  if (is_scotopic) {
    Fe = pow(10, (0.48 * kX));
    Fsc = mymin(1, (1 - pow(Pst / 124.4, 4)) / (1 - pow((OpticDia / OpticMag / 124.4), 4)));
    Fci = pow(10, (-0.4 * (1 - CIi / 2.0)));
    Fcb = pow(10, (-0.4 * (1 - CIb / 2.0)));
  } else {
     Fe = pow(10, (0.4 * kX));
     Fsc = mymin(1, pow((OpticDia / OpticMag / Pst), 2) * (1 - exp(-pow((Pst / 6.2), 2))) / (1 - exp(-pow((OpticDia / OpticMag / 6.2), 2))));
     Fci = 1;
     Fcb = 1;
  }
  Ft = 1 / OpticTrans;
  Fp = mymax(1, pow((Pst / (OpticMag * PupilDia(Age, Bback))), 2));
  Fa = pow((Pst / OpticDia), 2);
  Fr = (1 + 0.03 * pow((OpticMag * ObjectSize / CVA(Bback, SNi, helflag)), 2)) / pow(SNi, 2);
  Fm = pow(OpticMag, 2);
#if SWEHEL_DEBUG
  fprintf(stderr, "Pst=%f\n", Pst);
  fprintf(stderr, "Fb =%f\n", Fb);
  fprintf(stderr, "Fe =%f\n", Fe);
  fprintf(stderr, "Ft =%f\n", Ft);
  fprintf(stderr, "Fp =%f\n", Fp);
  fprintf(stderr, "Fa =%f\n", Fa);
  fprintf(stderr, "Fm =%f\n", Fm);
  fprintf(stderr, "Fsc=%f\n", Fsc);
  fprintf(stderr, "Fci=%f\n", Fci);
  fprintf(stderr, "Fcb=%f\n", Fcb);
  fprintf(stderr, "Fr =%f\n", Fr );
#endif
  if (TypeFactor == 0)
    return Fb * Fe * Ft * Fp * Fa * Fr * Fsc * Fci;
  else
    return Fb * Ft * Fp * Fa * Fm * Fsc * Fcb;
}

static int32 DeterObject(char *ObjectName)
{
  char s[AS_MAXCH];
  char *sp;
  int32 ipl;
  strcpy(s, ObjectName);
  for (sp = s; *sp != '\0'; sp++)
    *sp = tolower(*sp);
  if (strncmp(s, "sun", 3) == 0)
    return SE_SUN;
  if (strncmp(s, "venus", 5) == 0)
    return SE_VENUS;
  if (strncmp(s, "mars", 4) == 0)
    return SE_MARS;
  if (strncmp(s, "mercur", 6) == 0)
    return SE_MERCURY;
  if (strncmp(s, "jupiter", 7) == 0)
    return SE_JUPITER;
  if (strncmp(s, "saturn", 6) == 0)
    return SE_SATURN;
  if (strncmp(s, "uranus", 6) == 0)
    return SE_URANUS;
  if (strncmp(s, "neptun", 6) == 0)
    return SE_NEPTUNE;
  if (strncmp(s, "moon", 4) == 0)
    return SE_MOON;
  if ((ipl = atoi(s)) > 0) {
    ipl += SE_AST_OFFSET;
    return ipl;
  }
  return -1;
}

#if 0
int32 call_swe_calc(double tjd, int32 ipl, int32 iflag, double *x, char *serr)
{
  int32 retval = OK, ipli, i;
  double dtjd;
  static TLS double tjdsv[3];
  static TLS double xsv[3][6];
  static TLS int32 iflagsv[3];
  ipli = ipl;
  if (ipli > SE_MOON)
    ipli = 2;
  dtjd = tjd - tjdsv[ipli];
  if (tjdsv[ipli] != 0 && iflag == iflagsv[ipli] && fabs(dtjd) < 5.0 / 1440.0) {
    for (i = 0; i < 3; i++)
      x[i] = xsv[ipli][i] + dtjd * xsv[ipli][i+3];
    for (i = 3; i < 6; i++)
      x[i] = xsv[ipli][i];
  } else {
    retval = swe_calc(tjd, ipl, iflag, x, serr);
    tjdsv[ipli] = tjd;
    iflagsv[ipli] = iflag;
    for (i = 0; i < 6; i++)
      xsv[ipli][i] = x[i];
  }
  return retval;
}
#endif

static int32 call_swe_fixstar(char *star, double tjd, int32 iflag, double *xx, char *serr)
{
  int32 retval;
  char star2[AS_MAXCH];
  strcpy(star2, star);
  retval =  swe_fixstar(star2, tjd, iflag, xx, serr);
  return retval;
}

static int32 call_swe_fixstar_mag(char *star, double *mag, char *serr)
{
  int32 retval;
  char star2[AS_MAXCH];
  static TLS double dmag;
  static TLS char star_save[AS_MAXCH];
  if (strcmp(star, star_save) == 0) {
    *mag = dmag;
    return OK;
  }
  strcpy(star_save, star);
  strcpy(star2, star);
  retval = swe_fixstar_mag(star2, &dmag, serr);
  *mag = dmag;
  return retval;
}

static int32 call_swe_rise_trans(double tjd, int32 ipl, char *star, int32 helflag, int32 eventtype, double *dgeo, double atpress, double attemp, double *tret, char *serr)
{
  int32 retval;
  int32 iflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  char star2[AS_MAXCH];
  strcpy(star2, star);
  retval = swe_rise_trans(tjd, ipl, star2, iflag, eventtype, dgeo, atpress, attemp, tret, serr);
  return retval;
}

static int32 calc_rise_and_set(double tjd_start, int32 ipl, double *dgeo, double *datm, int32 eventflag, int32 helflag, double *trise, char *serr)
{
  int retc = OK, i;
  double sda, xs[6], xx[6], xaz[6], xaz2[6], dfac = 1/365.25;
  double rdi, rh;
  double tjd0 = tjd_start, tjdrise;
  double tjdnoon = (int) tjd0 - dgeo[0] / 15.0 / 24.0;
  int32 iflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  int32 epheflag = iflag;
  iflag |= SEFLG_EQUATORIAL;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT|SEFLG_TRUEPOS;
  if (swe_calc_ut(tjd0, SE_SUN, iflag, xs, serr) == 0) {
    if (serr != NULL)
      strcpy(serr, "error in calc_rise_and_set(): calc(sun) failed ");
    return ERR;
  }
  if (swe_calc_ut(tjd0, ipl, iflag, xx, serr) == 0) {
    if (serr != NULL)
      strcpy(serr, "error in calc_rise_and_set(): calc(sun) failed ");
    return ERR;
  }
  tjdnoon -= swe_degnorm(xs[0] - xx[0])/360.0 + 0;

  swe_azalt(tjd0, SE_EQU2HOR, dgeo, datm[0], datm[1], xx, xaz);
if (eventflag & SE_CALC_RISE) {
  if (xaz[2] > 0) {
    while (tjdnoon - tjd0 < 0.5) {tjdnoon += 1;}
    while (tjdnoon - tjd0 > 1.5) {tjdnoon -= 1;}
  } else {
    while (tjdnoon - tjd0 < 0.0) {tjdnoon += 1;}
    while (tjdnoon - tjd0 > 1.0) {tjdnoon -= 1;}
  }
} else {
  if (xaz[2] > 0) {
    while (tjd0 - tjdnoon > 0.5) { tjdnoon += 1;}
    while (tjd0 - tjdnoon < -0.5) { tjdnoon -= 1;}
  } else {
    while (tjd0 - tjdnoon > 0.0) { tjdnoon += 1;}
    while (tjd0 - tjdnoon < -1.0) { tjdnoon -= 1;}
  }
}

  if (swe_calc_ut(tjdnoon, ipl, iflag, xx, serr) == ERR) {
    if (serr != NULL)
      strcpy(serr, "error in calc_rise_and_set(): calc(sun) failed ");
    return ERR;
  }

  rdi = 0;
  if (ipl == SE_SUN)
    rdi = asin(696000000.0 / 1.49597870691e+11 / xx[2]) / DEGTORAD;
  else if (ipl == SE_MOON)
    rdi = asin(1737000.0 / 1.49597870691e+11 / xx[2]) / DEGTORAD;
  if (eventflag & SE_BIT_DISC_CENTER)
    rdi = 0;

  rh = -(34.5 / 60.0 + rdi);

  sda = acos(-tan(dgeo[1] * DEGTORAD) * tan(xx[1] * DEGTORAD)) * RADTODEG;

  if (eventflag & SE_CALC_RISE)
    tjdrise = tjdnoon - sda / 360.0;
  else
    tjdrise = tjdnoon + sda / 360.0;

  iflag = epheflag|SEFLG_SPEED|SEFLG_EQUATORIAL;
  if (ipl == SE_MOON)
    iflag |= SEFLG_TOPOCTR;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT|SEFLG_TRUEPOS;
  for (i = 0; i < 2; i++) {
    if (swe_calc_ut(tjdrise, ipl, iflag, xx, serr) == ERR) {

      return ERR;
    }
    swe_azalt(tjdrise, SE_EQU2HOR, dgeo, datm[0], datm[1], xx, xaz);
    xx[0] -= xx[3] * dfac;
    xx[1] -= xx[4] * dfac;
    swe_azalt(tjdrise - dfac, SE_EQU2HOR, dgeo, datm[0], datm[1], xx, xaz2);
    tjdrise -= (xaz[1] - rh) / (xaz[1] - xaz2[1]) *  dfac;

  }
  *trise = tjdrise;
  return retc;
}

static int32 my_rise_trans(double tjd, int32 ipl, char* starname, int32 eventtype, int32 helflag, double *dgeo, double *datm, double *tret, char *serr)
{
  int retc = OK;
  if (starname != NULL && *starname != '\0')
    ipl = DeterObject(starname);

  if (ipl != -1 && fabs(dgeo[1]) < 63) {
    retc = calc_rise_and_set(tjd, ipl, dgeo, datm, eventtype, helflag, tret, serr);

  } else {
    retc = call_swe_rise_trans(tjd, ipl, starname, helflag, eventtype, dgeo, datm[0], datm[1], tret, serr);
  }

  return retc;
}

static int32 RiseSet(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 RSEvent, int32 helflag, int32 Rim, double *tret, char *serr)
{
  int32 eventtype = RSEvent, Planet, retval;
  if (Rim == 0)
    eventtype |= SE_BIT_DISC_CENTER;
  Planet = DeterObject(ObjectName);
  if (Planet != -1)
    retval = my_rise_trans(JDNDaysUT, Planet, "", eventtype, helflag, dgeo, datm, tret, serr);
  else
    retval = my_rise_trans(JDNDaysUT, -1, ObjectName, eventtype, helflag, dgeo, datm, tret, serr);
  return retval;
}

static double SunRA(double JDNDaysUT, int32 helflag, char *serr)
{
  int imon, iday, iyar, calflag = SE_GREG_CAL;
  double dut;
  static TLS double tjdlast;
  static TLS double ralast;
  helflag += 0;
  *serr = '\0';
  if (JDNDaysUT == tjdlast)
    return ralast;
#ifndef SIMULATE_VICTORVB
  if (1) {
    double tjd_tt;
    double x[6];
    int32 epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
    int32 iflag = epheflag | SEFLG_EQUATORIAL;
    iflag |= SEFLG_NONUT | SEFLG_TRUEPOS;
    tjd_tt = JDNDaysUT + swe_deltat_ex(JDNDaysUT, epheflag, serr);
    if (swe_calc(tjd_tt, SE_SUN, iflag, x, serr) != ERR) {
      ralast = x[0];
      tjdlast = JDNDaysUT;
      return ralast;
    }
  }
#endif
  swe_revjul(JDNDaysUT, calflag, &iyar, &imon, &iday, &dut);
  tjdlast = JDNDaysUT;
  ralast = swe_degnorm((imon + (iday - 1) / 30.4 - 3.69) * 30);

  return ralast;
}

static double Kelvin(double Temp)
{

  return Temp + C2K;
}

static double TopoAltfromAppAlt(double AppAlt, double TempE, double PresE)
{
  double R = 0;
  double retalt = 0;
  if (AppAlt >= LowestAppAlt) {
    if (AppAlt > 17.904104638432)
      R = 0.97 / tan(AppAlt * DEGTORAD);
    else
      R = (34.46 + 4.23 * AppAlt + 0.004 * AppAlt * AppAlt) / (1 + 0.505 * AppAlt + 0.0845 * AppAlt * AppAlt);
    R = (PresE - 80) / 930 / (1 + 0.00008 * (R + 39) * (TempE - 10)) * R;
    retalt = AppAlt - R * Min2Deg;
  } else {
    retalt = AppAlt;
  }
  return retalt;
}

static double AppAltfromTopoAlt(double TopoAlt, double TempE, double PresE, int32 helflag)
{

  int i, nloop = 2;
  double newAppAlt = TopoAlt;
  double newTopoAlt = 0.0;
  double oudAppAlt = newAppAlt;
  double oudTopoAlt = newTopoAlt;
  double verschil, retalt;
  if (helflag & SE_HELFLAG_HIGH_PRECISION)
    nloop = 5;
  for (i = 0; i <= nloop; i++) {
    newTopoAlt = newAppAlt - TopoAltfromAppAlt(newAppAlt, TempE, PresE);

    verschil = newAppAlt - oudAppAlt;
    oudAppAlt = newTopoAlt - oudTopoAlt - verschil;
    if ((verschil != 0) && (oudAppAlt != 0))
      verschil = newAppAlt - verschil * (TopoAlt + newTopoAlt - newAppAlt) / oudAppAlt;
    else
      verschil = TopoAlt + newTopoAlt;
    oudAppAlt = newAppAlt;
    oudTopoAlt = newTopoAlt;
    newAppAlt = verschil;
  }
  retalt = TopoAlt + newTopoAlt;
  if (retalt < LowestAppAlt)
    retalt = TopoAlt;
  return retalt;
}

static double HourAngle(double TopoAlt, double TopoDecl, double Lat)
{
  double Alti = TopoAlt * DEGTORAD;
  double decli = TopoDecl * DEGTORAD;
  double Lati = Lat * DEGTORAD;
  double ha = (sin(Alti) - sin(Lati) * sin(decli)) / cos(Lati) / cos(decli);
  if (ha < -1) ha = -1;
  if (ha > 1) ha = 1;

  return acos(ha) / DEGTORAD / 15.0;
}

static int32 ObjectLoc(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 Angle, int32 helflag, double *dret, char *serr)
{
  double x[6], xin[3], xaz[3], tjd_tt;
  int32 Planet;
  int32 epheflag;
  int32 iflag = SEFLG_EQUATORIAL;
  epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  iflag |= epheflag;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT | SEFLG_TRUEPOS;
  if (Angle < 5) iflag = iflag | SEFLG_TOPOCTR;
  if (Angle == 7) Angle = 0;
  tjd_tt = JDNDaysUT + swe_deltat_ex(JDNDaysUT, epheflag, serr);
  Planet = DeterObject(ObjectName);
  if (Planet != -1) {
    if (swe_calc(tjd_tt, Planet, iflag, x, serr) == ERR)
      return ERR;
  } else {
    if (call_swe_fixstar(ObjectName, tjd_tt, iflag, x, serr) == ERR)
      return ERR;
  }
  if (Angle == 2 ||  Angle == 5) {
    *dret = x[1];
  } else {
    if (Angle == 3 || Angle == 6) {
      *dret = x[0];
    } else {
      xin[0] = x[0];
      xin[1] = x[1];
      swe_azalt(JDNDaysUT, SE_EQU2HOR, dgeo, datm[0], datm[1], xin, xaz);
      if (Angle == 0)
	*dret = xaz[1];
      if (Angle == 4)
	*dret = AppAltfromTopoAlt(xaz[1], datm[0], datm[1], helflag);
      if (Angle == 1) {
	xaz[0] += 180;
	if (xaz[0] >= 360)
	  xaz[0] -= 360;
        *dret = xaz[0];
      }
    }
  }
  return OK;
}

static int32 azalt_cart(double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 helflag, double *dret, char *serr)
{
  double x[6], xin[3], xaz[3], tjd_tt;
  int32 Planet;
  int32 epheflag;
  int32 iflag = SEFLG_EQUATORIAL;
  epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  iflag |= epheflag;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT | SEFLG_TRUEPOS;
  iflag = iflag | SEFLG_TOPOCTR;
  tjd_tt = JDNDaysUT + swe_deltat_ex(JDNDaysUT, epheflag, serr);
  Planet = DeterObject(ObjectName);
  if (Planet != -1) {
    if (swe_calc(tjd_tt, Planet, iflag, x, serr) == ERR)
      return ERR;
  } else {
    if (call_swe_fixstar(ObjectName, tjd_tt, iflag, x, serr) == ERR)
      return ERR;
  }
  xin[0] = x[0];
  xin[1] = x[1];
  swe_azalt(JDNDaysUT, SE_EQU2HOR, dgeo, datm[0], datm[1], xin, xaz);
  dret[0] = xaz[0];
  dret[1] = xaz[1];
  dret[2] = xaz[2];

  xaz[1] = xaz[2];
  xaz[2] = 1;
  swi_polcart(xaz, xaz);
  dret[3] = xaz[0];
  dret[4] = xaz[1];
  dret[5] = xaz[2];
  return OK;
}

static double DistanceAngle(double LatA, double LongA, double LatB, double LongB)
{
  double dlon = LongB - LongA;
  double dlat = LatB - LatA;

  double sindlat2 = sin(dlat / 2);
  double sindlon2 = sin(dlon / 2);
  double corde = sindlat2 * sindlat2 + cos(LatA) * cos(LatB) * sindlon2 *sindlon2;
  if (corde > 1) corde = 1;
  return 2 * asin(sqrt(corde));
}

static double kW(double HeightEye, double TempS, double RH)
{

  double WT = 0.031;
  WT *= 0.94 * (RH / 100.0) * exp(TempS / 15) * exp(-1 * HeightEye / scaleHwater);
  return WT;
}

static double kOZ(double AltS, double sunra, double Lat)
{
  double CHANGEKO, OZ, LT, kOZret;
  static TLS double koz_last, alts_last, sunra_last;
  double altslim = 0;
  if (AltS == alts_last && sunra == sunra_last)
    return koz_last;
  alts_last = AltS; sunra_last = sunra;
  OZ = 0.031;
  LT = Lat * DEGTORAD;

  kOZret = OZ * (3.0 + 0.4 * (LT * cos(sunra * DEGTORAD) - cos(3 * LT))) / 3.0;

  altslim = -AltS - 12;
  if (altslim < 0)
    altslim = 0;
  CHANGEKO = (100 - 11.6 * mymin(6, altslim)) / 100;
if ((0)) {
  static int a = 0;
  if (a == 0)
    printf("bsk=%f %f\n", kOZret, AltS);
  a = 1;
}
  koz_last = kOZret * CHANGEKO;
  return koz_last;
}

static double kR(double AltS, double HeightEye)
{

  double CHANGEK, LAMBDA;
  double val = -AltS - 12;
  if (val < 0) val = 0;
  if (val > 6) val = 6;

  CHANGEK = (1 - 0.166667 * val );
  LAMBDA = 0.55 + (CHANGEK - 1) * 0.04;

  return 0.1066 * exp(-1 * HeightEye / scaleHrayleigh) * pow(LAMBDA / 0.55 , -4);
}

static int Sgn(double x)
{
  if (x < 0)
    return -1;
  return 1;
}

static double ka(double AltS, double sunra, double Lat, double HeightEye, double TempS, double RH, double VR, char *serr)
{
  double CHANGEKA, LAMBDA, BetaVr, Betaa, kaact;
  double SL = Sgn(Lat);

  static TLS double alts_last, sunra_last, ka_last;
  if (AltS == alts_last && sunra == sunra_last)
    return ka_last;
  alts_last = AltS; sunra_last = sunra;
  CHANGEKA = (1 - 0.166667 * mymin(6, mymax(-AltS - 12, 0)));
  LAMBDA = 0.55 + (CHANGEKA - 1) * 0.04;
  if (VR != 0) {
    if (VR >= 1) {

      BetaVr = 3.912 / VR;
      Betaa = BetaVr - (kW(HeightEye, TempS, RH) / scaleHwater + kR(AltS, HeightEye) / scaleHrayleigh) * 1000 * astr2tau;
      kaact = Betaa * scaleHaerosol / 1000 * tau2astr;
      if (kaact < 0) {
	if (serr != NULL)
	  strcpy(serr, "The provided Meteorological range is too long, when taking into acount other atmospheric parameters");

      }
    } else {
      kaact = VR - kW(HeightEye, TempS, RH) - kR(AltS, HeightEye) - kOZ(AltS, sunra, Lat);
      if (kaact < 0) {
	if (serr != NULL)
	  strcpy(serr, "The provided atmosphic coeefficent (ktot) is too low, when taking into acount other atmospheric parameters");

      }
    }
  } else {

#ifdef SIMULATE_VICTORVB
    if (RH <= 0.00000001) RH = 0.00000001;
    if (RH >= 99.99999999) RH = 99.99999999;
#endif
    kaact = 0.1 * exp(-1 * HeightEye / scaleHaerosol) * pow(1 - 0.32 / log(RH / 100.0), 1.33) * (1 + 0.33 * SL * sin(sunra * DEGTORAD));
    kaact = kaact * pow(LAMBDA / 0.55, -1.3);
  }
  ka_last = kaact;
  return kaact;
}

static double kt(double AltS, double sunra, double Lat, double HeightEye, double TempS, double RH, double VR, int32 ExtType, char *serr)
{
  double kRact = 0;
  double kWact = 0;
  double kOZact = 0;
  double kaact = 0;
  if (ExtType == 2 || ExtType == 4)
    kRact = kR(AltS, HeightEye);
  if (ExtType == 1 || ExtType == 4)
    kWact = kW(HeightEye, TempS, RH);
  if (ExtType == 3 || ExtType == 4)
    kOZact = kOZ(AltS, sunra, Lat);
  if (ExtType == 0 || ExtType == 4)
    kaact = ka(AltS, sunra, Lat, HeightEye, TempS, RH, VR, serr);
  if (kaact < 0)
    kaact = 0;
  return kWact + kRact + kOZact + kaact;
}

static double Airmass(double AppAltO, double Press)
{
  double airm, zend;
  zend = (90 - AppAltO) * DEGTORAD;
  if (zend > PI / 2)
    zend = PI / 2;
  airm = 1 / (cos(zend) + 0.025 * exp(-11 * cos(zend)));
  return Press / 1013 * airm;
}

static double Xext(double scaleH, double zend, double Press)
{
  return Press / 1013.0 / (cos(zend) + 0.01 * sqrt(scaleH / 1000.0) * exp(-30.0 / sqrt(scaleH / 1000.0) * cos(zend)));
}

static double Xlay(double scaleH, double zend, double Press)
{

  double a = sin(zend) / (1.0 + (scaleH / Ra));
  return Press / 1013.0 / sqrt(1.0 - a * a);
}

static double TempEfromTempS(double TempS, double HeightEye, double Lapse)
{
  return TempS - Lapse * HeightEye;
}

static double PresEfromPresS(double TempS, double Press, double HeightEye)
{
  return Press * exp(-9.80665 * 0.0289644 / (Kelvin(TempS) + 3.25 * HeightEye / 1000) / 8.31441 * HeightEye);
}

static double Deltam(double AltO, double AltS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)
{
  double zend, xR, XW, Xa, XOZ;
  double PresE = PresEfromPresS(datm[1], datm[0], HeightEye);
  double TempE = TempEfromTempS(datm[1], HeightEye, LapseSA);
  double AppAltO = AppAltfromTopoAlt(AltO, TempE, PresE, helflag);
  double deltam;
  static TLS double alts_last, alto_last, sunra_last, deltam_last;
  if (AltS == alts_last && AltO == alto_last && sunra == sunra_last)
    return deltam_last;
  alts_last = AltS; alto_last = AltO; sunra_last = sunra;
  if (staticAirmass == 0) {
    zend = (90 - AppAltO) * DEGTORAD;
    if (zend > PI / 2)
      zend = PI / 2;

    xR = Xext(scaleHrayleigh, zend, datm[0]);
    XW = Xext(scaleHwater, zend, datm[0]);
    Xa = Xext(scaleHaerosol, zend, datm[0]);
    XOZ = Xlay(scaleHozone, zend, datm[0]);
    deltam = kR(AltS, HeightEye) * xR + kt(AltS, sunra, Lat, HeightEye, datm[1], datm[2], datm[3], 0, serr) * Xa + kOZ(AltS, sunra, Lat) * XOZ + kW(HeightEye, datm[1], datm[2]) * XW;
  } else {
    deltam = kt(AltS, sunra, Lat, HeightEye, datm[1], datm[2], datm[3], 4, serr) * Airmass(AppAltO, datm[0]);
  }
  deltam_last = deltam;
  return deltam;
}

static double Bn(double AltO, double JDNDayUT, double AltS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)
{
  double PresE = PresEfromPresS(datm[1], datm[0], HeightEye);
  double TempE = TempEfromTempS(datm[1], HeightEye, LapseSA);
  double AppAltO = AppAltfromTopoAlt(AltO, TempE, PresE, helflag);
  double zend, YearB, MonthB, DayB, Bna, kX, Bnb;
  double B0 = 0.0000000000001, dut;
  int iyar, imon, iday;

  if (AppAltO < 10)
    AppAltO = 10;
  zend = (90 - AppAltO) * DEGTORAD;

  swe_revjul(JDNDayUT, SE_GREG_CAL, &iyar, &imon, &iday, &dut);
  YearB = iyar; MonthB = imon; DayB = iday;
  Bna = B0 * (1 + 0.3 * cos(6.283 * (YearB + ((DayB - 1) / 30.4 + MonthB - 1) / 12 - 1990.33) / 11.1));
  kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);

  Bnb = Bna * (0.4 + 0.6 / sqrt(1 - 0.96 * pow(sin(zend), 2))) * pow(10, -0.4 * kX);
  return mymax(Bnb, 0) * erg2nL;
}

static int32 Magnitude(double JDNDaysUT, double *dgeo, char *ObjectName, int32 helflag, double *dmag, char *serr)
{
  double x[20];
  int32 Planet, iflag, epheflag;
  epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  *dmag = -99.0;
  Planet = DeterObject(ObjectName);
  iflag = SEFLG_TOPOCTR | SEFLG_EQUATORIAL | epheflag;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT|SEFLG_TRUEPOS;
  if (Planet != -1) {

    swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);
    if (swe_pheno_ut(JDNDaysUT, Planet, iflag, x, serr) == ERR)
      return ERR;
    *dmag = x[4];
  } else {
    if (call_swe_fixstar_mag(ObjectName, dmag, serr) == ERR)
      return ERR;
  }
  return OK;
}

#if 0
static int32 fast_magnitude(double tjd, double *dgeo, char *ObjectName, int32 helflag, double *dmag, char *serr)
{
  int32 retval = OK, ipl, ipli;
  double dtjd;
  static TLS double tjdsv[3];
  static TLS double dmagsv[3];
  static TLS int32 helflagsv[3];
  ipl = DeterObject(ObjectName);
  ipli = ipl;
  if (ipli > SE_MOON)
    ipli = 2;
  dtjd = tjd - tjdsv[ipli];
  if (tjdsv[ipli] != 0 && helflag == helflagsv[ipli] && fabs(dtjd) < 5.0 / 1440.0) {
    *dmag = dmagsv[ipli];
  } else {
    retval = Magnitude(tjd, dgeo, ObjectName, helflag, dmag, serr);
    tjdsv[ipli] = tjd;
    helflagsv[ipli] = helflag;
    dmagsv[ipli] = *dmag;
  }
  return retval;
}
#endif

static double MoonsBrightness(double dist, double phasemoon)
{
  double log10 = 2.302585092994;

  return -21.62 + 5 * log(dist / (Ra / 1000)) / log10 + 0.026 * fabs(phasemoon) + 0.000000004 * pow(phasemoon, 4);
}

static double MoonPhase(double AltM, double AziM, double AltS, double AziS)
{
  double AltMi = AltM * DEGTORAD;
  double AltSi = AltS * DEGTORAD;
  double AziMi = AziM * DEGTORAD;
  double AziSi = AziS * DEGTORAD;
  double MoonAvgPar = 0.95;

  return 180 - acos(cos(AziSi - AziMi - MoonAvgPar * DEGTORAD) * cos(AltMi + MoonAvgPar * DEGTORAD) * cos(AltSi) + sin(AltSi) * sin(AltMi + MoonAvgPar * DEGTORAD)) / DEGTORAD;
}

static double Bm(double AltO, double AziO, double AltM, double AziM, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)
{
  double M0 = -11.05;
  double Bm = 0;
  double RM, kXM, kX, C3, FM, phasemoon, MM;
  double lunar_radius = 0.25 * DEGTORAD;
  AS_BOOL object_is_moon = FALSE;
  if (AltO == AltM && AziO == AziM)
    object_is_moon = TRUE;
  if (AltM > -0.26 && !object_is_moon) {

    RM = DistanceAngle(AltO * DEGTORAD, AziO * DEGTORAD, AltM * DEGTORAD, AziM * DEGTORAD) / DEGTORAD;
    if (RM <= lunar_radius)
      RM = lunar_radius;
    kXM = Deltam(AltM, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
    kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
    C3 = pow(10, -0.4 * kXM);
    FM = (62000000.0) / RM / RM + pow(10, 6.15 - RM / 40) + pow(10, 5.36) * (1.06 + pow(cos(RM * DEGTORAD), 2));
    Bm = FM * C3 + 440000 * (1 - C3);
    phasemoon = MoonPhase(AltM, AziM, AltS, AziS);
    MM = MoonsBrightness(MoonDistance, phasemoon);
    Bm = Bm * pow(10, -0.4 * (MM - M0 + 43.27));
    Bm = Bm * (1 - pow(10, -0.4 * kX));
  }
  Bm = mymax(Bm, 0) * erg2nL;
  return Bm;
}

static double Btwi(double AltO, double AziO, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)
{
double M0 = -11.05;
double MS = -26.74;
double PresE = PresEfromPresS(datm[1], datm[0], HeightEye);
double TempE = TempEfromTempS(datm[1], HeightEye, LapseSA);
double AppAltO = AppAltfromTopoAlt(AltO, TempE, PresE, helflag);
double ZendO = 90 - AppAltO;
double RS = DistanceAngle(AltO * DEGTORAD, AziO * DEGTORAD, AltS * DEGTORAD, AziS * DEGTORAD) / DEGTORAD;
double kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
double k = kt(AltS, sunra, Lat, HeightEye, datm[1], datm[2], datm[3], 4, serr);

double Btwi = pow(10, -0.4 * (MS - M0 + 32.5 - AltS - (ZendO / (360 * k))));
Btwi = Btwi * (100 / RS) * (1 - pow(10, -0.4 * kX));
Btwi = mymax(Btwi, 0) * erg2nL;
return Btwi;
}

static double Bday(double AltO, double AziO, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)
{
  double M0 = -11.05;
  double MS = -26.74;
  double RS = DistanceAngle(AltO * DEGTORAD, AziO * DEGTORAD, AltS * DEGTORAD, AziS * DEGTORAD) / DEGTORAD;
  double kXS = Deltam(AltS, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
  double kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);

  double C4 = pow(10, -0.4 * kXS);
  double FS = (62000000.0) / RS / RS + pow(10, (6.15 - RS / 40)) + pow(10, 5.36) * (1.06 + pow(cos(RS * DEGTORAD), 2));
  double Bday = FS * C4 + 440000.0 * (1 - C4);
  Bday = Bday * pow(10, (-0.4 * (MS - M0 + 43.27)));
  Bday = Bday * (1 - pow(10, -0.4 * kX));
  Bday = mymax(Bday, 0) * erg2nL;
  return Bday;
}

static double Bcity(double Value, double Press)
{
  double Bcity = Value;
  Press += 0.0;
  Bcity = mymax(Bcity, 0);
  return Bcity;
}

static double Bsky(double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, char *serr)
{
  double Bsky = 0;
  if (AltS < -3) {
    Bsky += Btwi(AltO, AziO, AltS, AziS, sunra, Lat, HeightEye, datm, helflag, serr);
  } else {
    if (AltS > 4) {
      Bsky += Bday(AltO, AziO, AltS, AziS, sunra, Lat, HeightEye, datm, helflag, serr);
    } else {
      Bsky += mymin(Bday(AltO, AziO, AltS, AziS, sunra, Lat, HeightEye, datm, helflag, serr), Btwi(AltO, AziO, AltS, AziS, sunra, Lat, HeightEye, datm, helflag, serr));
if ((0)) {
  static int a = 0;
  if (a == 0)
    printf("bsk=%f\n", Bsky);
  a = 1;
}
    }
  }

  if (Bsky < 200000000.0)
    Bsky += Bm(AltO, AziO, AltM, AziM, AltS, AziS, sunra, Lat, HeightEye, datm, helflag, serr);
  if (AltS <= 0)
    Bsky += Bcity(0, datm[0]);

  if (Bsky < 5000)
    Bsky = Bsky + Bn(AltO, JDNDaysUT, AltS, sunra, Lat, HeightEye, datm, helflag, serr);

  return Bsky;
}

static void default_heliacal_parameters(double *datm, double *dgeo, double *dobs, int helflag)
{
  int i;
  if (datm[0] <= 0) {

    datm[0] = 1013.25 * pow(1 - 0.0065 * dgeo[2] / 288, 5.255);

    if (datm[1] == 0)
      datm[1] = 15 - 0.0065 * dgeo[2];

    if (datm[2] == 0)
      datm[2] = 40;

  } else {
#ifndef SIMULATE_VICTORVB
    if (datm[2] <= 0.00000001) datm[2] = 0.00000001;
    if (datm[2] >= 99.99999999) datm[2] = 99.99999999;
#endif
  }

  if (dobs[0] == 0)
    dobs[0] = 36;

  if (dobs[1] == 0)
    dobs[1] = 1;
  if (!(helflag & SE_HELFLAG_OPTICAL_PARAMS)) {
    for (i = 2; i <= 5; i++)
      dobs[i] = 0;
  }

  if (dobs[3] == 0) {
    dobs[2] = 1;
    dobs[3] = 1;

  }
}

static double VisLimMagn(double *dobs, double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AltS, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, int32 *scotopic_flag, char *serr)
{
  double C1, C2, Th, kX, Bsk, CorrFactor1, CorrFactor2;
  double log10 = 2.302585092994;
  AS_BOOL is_scotopic = FALSE;

  Bsk = Bsky(AltO, AziO, AltM, AziM, JDNDaysUT, AltS, AziS, sunra, Lat, HeightEye, datm, helflag, serr);

  kX = Deltam(AltO, AltS, sunra, Lat, HeightEye, datm, helflag, serr);
if ((0)) {
  static int a = 0;
  if (a == 0)
    printf("bsk=%f, kx=%f\n", Bsk, kX);
  a = 1;
}

  CorrFactor1 = OpticFactor(Bsk, kX, dobs, JDNDaysUT, "", 1, helflag);
  CorrFactor2 = OpticFactor(Bsk, kX, dobs, JDNDaysUT, "", 0, helflag);

  if (Bsk < 1645)
    is_scotopic = TRUE;
  if (helflag & SE_HELFLAG_VISLIM_PHOTOPIC)
    is_scotopic = FALSE;
  if (helflag & SE_HELFLAG_VISLIM_SCOTOPIC)
    is_scotopic = TRUE;

  if (is_scotopic) {
    C1 = 1.5848931924611e-10;
    C2 = 0.012589254117942;
    if (scotopic_flag != NULL)
      *scotopic_flag = 1;
  } else {
    C1 = 4.4668359215096e-9;
    C2 = 1.2589254117942e-6;
    if (scotopic_flag != NULL)
      *scotopic_flag = 0;
  }
  if (scotopic_flag != NULL) {
    if (BNIGHT * BNIGHT_FACTOR > Bsk && BNIGHT / BNIGHT_FACTOR < Bsk)
      *scotopic_flag |= 2;
  }

  Bsk = Bsk * CorrFactor1;
  Th = C1 * pow(1 + sqrt(C2 * Bsk), 2) * CorrFactor2;
#if SWEHEL_DEBUG
  fprintf(stderr, "Bsk=%f, ", Bsk);
  fprintf(stderr, "kX =%f, ", kX);
  fprintf(stderr, "Th =%f, ", Th);
  fprintf(stderr, "CorrFactor1=%f, ", CorrFactor1);
  fprintf(stderr, "CorrFactor2=%f\n", CorrFactor2);
#endif

#if 0
  if (SN <= 0.00000001)
    SN = 0.00000001;
  return -16.57 - 2.5 * (log(Th) / log10) - kX + 5.0 * (log(SN) / log10);
#endif
  return -16.57 - 2.5 * (log(Th) / log10);
}

static char *tolower_string_star(char *str)
{
  char *sp;
  for (sp = str; *sp != '\0' && *sp != ','; sp++)
    *sp = tolower(*sp);
  return str;
}

int32 CALL_CONV swe_vis_limit_mag(double tjdut, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 helflag, double *dret, char *serr)
{
  int32 retval = OK, i, scotopic_flag = 0;
  double AltO, AziO, AltM, AziM, AltS, AziS;
  double sunra;
  for (i = 0; i < 7; i++)
    dret[i] = 0;
  tolower_string_star(ObjectName);
  if (DeterObject(ObjectName) == SE_SUN) {
    if (serr != NULL) {
      strcpy(serr, "it makes no sense to call swe_vis_limit_mag() for the Sun");
    }
    return ERR;
  }
  swi_set_tid_acc(tjdut, helflag, 0, serr);
  sunra = SunRA(tjdut, helflag, serr);
  default_heliacal_parameters(datm, dgeo, dobs, helflag);
  swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);
  if (ObjectLoc(tjdut, dgeo, datm, ObjectName, 0, helflag, &AltO, serr) == ERR)
    return ERR;
  if (AltO < 0) {
    if (serr != NULL)
      strcpy(serr, "object is below local horizon");
    *dret = -100;
    return -2;
  }
  if (ObjectLoc(tjdut, dgeo, datm, ObjectName, 1, helflag, &AziO, serr) == ERR)
    return ERR;
  if (helflag & SE_HELFLAG_VISLIM_DARK) {
    AltS = -90;
    AziS = 0;
  } else {
    if (ObjectLoc(tjdut, dgeo, datm, "sun", 0, helflag, &AltS, serr) == ERR)
      return ERR;
    if (ObjectLoc(tjdut, dgeo, datm, "sun", 1, helflag, &AziS, serr) == ERR)
      return ERR;
  }
  if (strncmp(ObjectName, "moon", 4) == 0 ||
      (helflag & SE_HELFLAG_VISLIM_DARK) ||
      (helflag & SE_HELFLAG_VISLIM_NOMOON)
     ) {
    AltM = -90; AziM = 0;
  } else {
    if (ObjectLoc(tjdut, dgeo, datm, "moon", 0, helflag, &AltM, serr) == ERR)
      return ERR;
    if (ObjectLoc(tjdut, dgeo, datm, "moon", 1, helflag, &AziM, serr) == ERR)
      return ERR;
  }
#if SWEHEL_DEBUG
{
  int i;
  for (i = 0; i < 6;i++)
    printf("dobs[%d] = %f\n", i, dobs[i]);
  printf("AltO = %.10f, AziO = %.10f\n", AltO, AziO);
  printf("AltM = %.10f, AziM = %.10f\n", AltM, AziM);
  printf("AltS = %.10f, AziS = %.10f\n", AltS, AziS);
  printf("JD = %.10f\n", tjdut);
  printf("lat = %f, eyeh = %f\n", dgeo[1], dgeo[2]);
  for (i = 0; i < 4;i++)
    printf("datm[%d] = %f\n", i, datm[i]);
  printf("helflag = %d\n", helflag);
}
#endif
  dret[0] = VisLimMagn(dobs, AltO, AziO, AltM, AziM, tjdut, AltS, AziS, sunra, dgeo[1], dgeo[2], datm, helflag, &scotopic_flag, serr);
  dret[1] = AltO;
  dret[2] = AziO;
  dret[3] = AltS;
  dret[4] = AziS;
  dret[5] = AltM;
  dret[6] = AziM;
  if (Magnitude(tjdut, dgeo, ObjectName, helflag, &(dret[7]), serr) == ERR)
    return ERR;
  retval = scotopic_flag;

  return retval;
}

static int32 TopoArcVisionis(double Magn, double *dobs, double AltO, double AziO, double AltM, double AziM, double JDNDaysUT, double AziS, double sunra, double Lat, double HeightEye, double *datm, int32 helflag, double *dret, char *serr)
{
  double Xm, Ym, AltSi, AziSi;
  double xR = 0;
  double Xl = 45;
  double Yl, Yr;
  Yl = Magn - VisLimMagn(dobs, AltO, AziO, AltM, AziM, JDNDaysUT, AltO - Xl, AziS, sunra, Lat, HeightEye, datm, helflag, NULL, serr);

  Yr = Magn - VisLimMagn(dobs, AltO, AziO, AltM, AziM, JDNDaysUT, AltO - xR, AziS, sunra, Lat, HeightEye, datm, helflag, NULL, serr);

  if ((Yl * Yr) <= 0) {
    while(fabs(xR - Xl) > epsilon) {

      Xm = (xR + Xl) / 2.0;
      AltSi = AltO - Xm;
      AziSi = AziS;
      Ym = Magn - VisLimMagn(dobs, AltO, AziO, AltM, AziM, JDNDaysUT, AltSi, AziSi, sunra, Lat, HeightEye, datm, helflag, NULL, serr);

      if ((Yl * Ym) > 0) {

	Xl = Xm;
	Yl = Ym;
      } else {

	xR = Xm;
	Yr = Ym;
      }
    }
    Xm = (xR + Xl) / 2.0;
  } else {
    Xm = 99;
  }
  if (Xm < AltO)
    Xm = AltO;
  *dret = Xm;
  return OK;
}

int32 CALL_CONV swe_topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr)
{
  double sunra;
  swi_set_tid_acc(tjdut, helflag, 0, serr);
  sunra = SunRA(tjdut, helflag, serr);
  if (serr != NULL && *serr != '\0')
    return ERR;
  default_heliacal_parameters(datm, dgeo, dobs, helflag);
  return TopoArcVisionis(mag, dobs, alt_obj, azi_obj, alt_moon, azi_moon, tjdut, azi_sun, sunra, dgeo[1], dgeo[2], datm, helflag, dret, serr);
}

static int32 HeliacalAngle(double Magn, double *dobs, double AziO, double AltM, double AziM, double JDNDaysUT, double AziS, double *dgeo, double *datm, int32 helflag, double *dangret, char *serr)
{
  double x, minx, maxx, xmin, ymin, Xl, xR, Yr, Yl, Xm, Ym, xmd, ymd;
  double Arc, DELTAx;
  double sunra = SunRA(JDNDaysUT, helflag, serr);
  double Lat = dgeo[1];
  double HeightEye = dgeo[2];
  if (PLSV == 1) {
    dangret[0] = criticalangle;
    dangret[1] = criticalangle + Magn * 2.492 + 13.447;
    dangret[2] = -(Magn * 2.492 + 13.447);
    return OK;
  }
  minx = 2;
  maxx = 20;
  xmin = 0;
  ymin = 10000;
  for (x = minx; x <= maxx; x++) {
    if (TopoArcVisionis(Magn, dobs, x, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, Lat, HeightEye, datm, helflag, &Arc, serr) == ERR)
      return ERR;
    if (Arc < ymin) {
      ymin = Arc;
      xmin = x;
    }
  }
  Xl = xmin - 1;
  xR = xmin + 1;
  if (TopoArcVisionis(Magn, dobs, xR, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, Lat, HeightEye, datm, helflag, &Yr, serr) == ERR)
    return ERR;
  if (TopoArcVisionis(Magn, dobs, Xl, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, Lat, HeightEye, datm, helflag, &Yl, serr) == ERR)
    return ERR;

  while(fabs(xR - Xl) > 0.1) {

    Xm = (xR + Xl) / 2.0;
    DELTAx = 0.025;
    xmd = Xm + DELTAx;
    if (TopoArcVisionis(Magn, dobs, Xm, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, Lat, HeightEye, datm, helflag, &Ym, serr) == ERR)
      return ERR;
    if (TopoArcVisionis(Magn, dobs, xmd, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, Lat, HeightEye, datm, helflag, &ymd, serr) == ERR)
      return ERR;
    if (Ym >= ymd) {

      Xl = Xm;
      Yl = Ym;
    } else {

      xR = Xm;
      Yr = Ym;
    }
  }
  Xm = (xR + Xl) / 2.0;
  Ym = (Yr + Yl) / 2.0;
  dangret[1] = Ym;
  dangret[2] = Xm - Ym;
  dangret[0] = Xm;
  return OK;
}

int32 CALL_CONV swe_heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr)
{
  if (dgeo[2] < SEI_ECL_GEOALT_MIN || dgeo[2] > SEI_ECL_GEOALT_MAX) {
    if (serr != NULL)
      sprintf(serr, "location for heliacal events must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  swi_set_tid_acc(tjdut, helflag, 0, serr);
  default_heliacal_parameters(datm, dgeo, dobs, helflag);
  return HeliacalAngle(mag, dobs, azi_obj, alt_moon, azi_moon, tjdut, azi_sun, dgeo, datm, helflag, dret, serr);
}

static double WidthMoon(double AltO, double AziO, double AltS, double AziS, double parallax)
{

  double GeoAltO = AltO + parallax;
  return 0.27245 * parallax * (1 + sin(GeoAltO * DEGTORAD) * sin(parallax * DEGTORAD)) * (1 - cos((AltS - GeoAltO) * DEGTORAD) * cos((AziS - AziO) * DEGTORAD));
}

static double LengthMoon(double W, double Diamoon)
{
  double Wi, D;
  if (Diamoon == 0) Diamoon = AvgRadiusMoon * 2;
  Wi = W * 60;
  D = Diamoon * 60;

  return (D - 0.3 * (D + Wi) / 2.0 / Wi) / 60.0;
}

static double qYallop(double W, double GeoARCVact)
{
  double Wi = W * 60;
  return (GeoARCVact - (11.8371 - 6.3226 * Wi + 0.7319 * Wi * Wi - 0.1018 * Wi * Wi * Wi)) / 10;
}

static double crossing(double A, double B, double C, double D)
{
return (C - A) / ((B - A) - (D - C));
}

static int32 DeterTAV(double *dobs, double JDNDaysUT, double *dgeo, double *datm, char *ObjectName, int32 helflag, double *dret, char *serr)
{
  double Magn, AltO, AziS, AziO, AziM, AltM;
  double sunra = SunRA(JDNDaysUT, helflag, serr);
  if (Magnitude(JDNDaysUT, dgeo, ObjectName, helflag, &Magn, serr) == ERR)
    return ERR;
  if (ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 0, helflag, &AltO, serr) == ERR)
    return ERR;
  if (ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 1, helflag, &AziO, serr) == ERR)
    return ERR;
  if (strncmp(ObjectName, "moon", 4) == 0) {
    AltM = -90;
    AziM = 0;
  } else {
    if (ObjectLoc(JDNDaysUT, dgeo, datm, "moon", 0, helflag, &AltM, serr) == ERR)
      return ERR;
    if (ObjectLoc(JDNDaysUT, dgeo, datm, "moon", 1, helflag, &AziM, serr) == ERR)
      return ERR;
  }
  if (ObjectLoc(JDNDaysUT, dgeo, datm, "sun", 1, helflag, &AziS, serr) == ERR)
    return ERR;
  if (TopoArcVisionis(Magn, dobs, AltO, AziO, AltM, AziM, JDNDaysUT, AziS, sunra, dgeo[1], dgeo[2], datm, helflag, dret, serr) == ERR)
    return ERR;
  return OK;
}

static double x2min(double A, double B, double C)
{
  double term = A + C - 2 * B;
  if (term == 0)
    return 0;
  return -(A - C) / 2.0 / term;
}

static double funct2(double A, double B, double C, double x)
{
  return (A + C - 2 * B) / 2.0 * x * x + (A - C) / 2.0 * x + B;
}

static void strcpy_VBsafe(char *sout, char *sin)
{
  char *sp, *sp2;
  int iw = 0;
  sp = sin;
  sp2 = sout;

  while((isalnum((int) *sp) || *sp == ' ' || *sp == '-' || *sp == ',') && iw < 30) {
    *sp2 = *sp;
    sp++; sp2++; iw++;
  }
  *sp2 = '\0';
}

int32 CALL_CONV swe_heliacal_pheno_ut(double JDNDaysUT, double *dgeo, double *datm, double *dobs, char *ObjectNameIn, int32 TypeEvent, int32 helflag, double *darr, char *serr)
{
  double AziS = 0, AltS = 0, AltS2 = 0, AziO = 0, AltO = 0, AltO2 = 0;
  double GeoAltO = 0, AppAltO = 0, DAZact = 0, TAVact = 0, ParO = 0, MagnO = 0;
  double ARCVact, ARCLact, kact, WMoon, LMoon = 0, qYal, qCrit;
  double RiseSetO, RiseSetS, Lag, TbYallop, TfirstVR, TlastVR, TbVR;
  double MinTAV = 0, MinTAVact, Ta, Tc, TimeStep, TimePointer, MinTAVoud = 0, DeltaAltoud = 0, DeltaAlt, TvisVR, crosspoint;
  double OldestMinTAV, extrax, illum;
  double elong, attr[30];
  double TimeCheck, LocalminCheck;
  int32 retval = OK, RS, Planet;
  AS_BOOL noriseO = FALSE;
  char ObjectName[AS_MAXCH];
  double sunra;
  int32 iflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  if (dgeo[2] < SEI_ECL_GEOALT_MIN || dgeo[2] > SEI_ECL_GEOALT_MAX) {
    if (serr != NULL)
      sprintf(serr, "location for heliacal events must be between %.0f and %.0f m above sea", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  swi_set_tid_acc(JDNDaysUT, helflag, 0, serr);
  sunra = SunRA(JDNDaysUT, helflag, serr);

  strcpy_VBsafe(ObjectName, ObjectNameIn);
  tolower_string_star(ObjectName);
  default_heliacal_parameters(datm, dgeo, dobs, helflag);
  swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);
  retval = ObjectLoc(JDNDaysUT, dgeo, datm, "sun", 1, helflag, &AziS, serr);
  if (retval == OK)
    retval = ObjectLoc(JDNDaysUT, dgeo, datm, "sun", 0, helflag, &AltS, serr);
  if (retval == OK)
    retval = ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 1, helflag, &AziO, serr);
  if (retval == OK)
    retval = ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 0, helflag, &AltO, serr);
  if (retval == OK)
    retval = ObjectLoc(JDNDaysUT, dgeo, datm, ObjectName, 7, helflag, &GeoAltO, serr);
  if (retval == ERR)
    return ERR;
  AppAltO = AppAltfromTopoAlt(AltO, datm[1], datm[0], helflag);
  DAZact = AziS - AziO;
  TAVact = AltO - AltS;

  ParO = GeoAltO - AltO;
  if (Magnitude(JDNDaysUT, dgeo, ObjectName, helflag, &MagnO, serr) == ERR)
    return ERR;
  ARCVact = TAVact + ParO;
  ARCLact = acos(cos(ARCVact * DEGTORAD) * cos(DAZact * DEGTORAD)) / DEGTORAD;
  Planet = DeterObject(ObjectName);
  if (Planet == -1) {
    elong = ARCLact;
    illum = 100;
  } else {
    retval = swe_pheno_ut(JDNDaysUT, Planet, iflag|(SEFLG_TOPOCTR|SEFLG_EQUATORIAL), attr, serr);
    if (retval == ERR) return ERR;
    elong = attr[2];
    illum = attr[1] * 100;
  }
  kact = kt(AltS, sunra, dgeo[1], dgeo[2], datm[1], datm[2], datm[3], 4, serr);
  if ((0)) {
darr[26] = kR(AltS, dgeo[2]);
darr[27] = kW(dgeo[2], datm[1], datm[2]);
darr[28] = kOZ(AltS, sunra, dgeo[1]);
darr[29] = ka(AltS, sunra, dgeo[1], dgeo[2], datm[1], datm[2], datm[3], serr);
darr[30] = darr[26] + darr[27] + darr[28] + darr[29];
  }
  WMoon = 0;
  qYal = 0;
  qCrit = 0;
  LMoon = 0;
  if (Planet == SE_MOON) {
    WMoon = WidthMoon(AltO, AziO, AltS, AziS, ParO);
    LMoon = LengthMoon(WMoon, 0);
    qYal = qYallop(WMoon, ARCVact);
    if (qYal > 0.216) qCrit = 1;
    if (qYal < 0.216 && qYal > -0.014) qCrit = 2;
    if (qYal < -0.014 && qYal > -0.16) qCrit = 3;
    if (qYal < -0.16 && qYal > -0.232) qCrit = 4;
    if (qYal < -0.232 && qYal > -0.293) qCrit = 5;
    if (qYal < -0.293) qCrit = 6;
  }

  RS = 2;
  if (TypeEvent == 1 || TypeEvent == 4) RS = 1;
  retval = RiseSet(JDNDaysUT - 4.0 / 24.0, dgeo, datm, "sun", RS, helflag, 0, &RiseSetS, serr);
  if (retval == ERR)
    return ERR;
  retval = RiseSet(JDNDaysUT - 4.0 / 24.0, dgeo, datm, ObjectName, RS, helflag, 0, &RiseSetO, serr);
  if (retval == ERR)
    return ERR;
  TbYallop = TJD_INVALID;
  if (retval == -2) {
    Lag = 0;
    noriseO = TRUE;
  } else {
    Lag = RiseSetO - RiseSetS;
    if (Planet == SE_MOON)
      TbYallop = (RiseSetO * 4 + RiseSetS * 5) / 9.0;
  }
  if ((TypeEvent == 3 || TypeEvent == 4) && (Planet == -1 || Planet >= SE_MARS)) {
    TfirstVR = TJD_INVALID;
    TbVR = TJD_INVALID;
    TlastVR = TJD_INVALID;
    TvisVR = 0;
    MinTAV = 0;
    goto output_heliacal_pheno;
  }

    MinTAVact = 199;
    DeltaAlt = 0;
    OldestMinTAV = 0;
    Ta = 0;
    Tc = 0;
    TbVR = 0;
    TimeStep = -TimeStepDefault / 24.0 / 60.0;
    if (RS == 2) TimeStep = -TimeStep;
    TimePointer = RiseSetS - TimeStep;
    do {
      TimePointer = TimePointer + TimeStep;
      OldestMinTAV = MinTAVoud;
      MinTAVoud = MinTAVact;
      DeltaAltoud = DeltaAlt;
      retval = ObjectLoc(TimePointer, dgeo, datm, "sun", 0, helflag, &AltS2, serr);
      if (retval == OK)
	retval = ObjectLoc(TimePointer, dgeo, datm, ObjectName, 0, helflag, &AltO2, serr);
      if (retval != OK)
        return ERR;
      DeltaAlt = AltO2 - AltS2;
      if (DeterTAV(dobs, TimePointer, dgeo, datm, ObjectName, helflag, &MinTAVact, serr) == ERR)
        return ERR;
      if (MinTAVoud < MinTAVact && TbVR == 0) {

	TimeCheck = TimePointer + Sgn(TimeStep) * LocalMinStep / 24.0 / 60.0;
	if (RiseSetO != 0) {
	  if (TimeStep > 0)
	    TimeCheck = mymin(TimeCheck, RiseSetO);
	  else
	    TimeCheck = mymax(TimeCheck, RiseSetO);
	}
        if (DeterTAV(dobs, TimeCheck, dgeo, datm, ObjectName, helflag, &LocalminCheck, serr) == ERR)
	  return ERR;
	if (LocalminCheck > MinTAVact) {
	  extrax = x2min(MinTAVact, MinTAVoud, OldestMinTAV);
	  TbVR = TimePointer - (1 - extrax) * TimeStep;
	  MinTAV = funct2(MinTAVact, MinTAVoud, OldestMinTAV, extrax);
	}
      }
      if (DeltaAlt > MinTAVact && Tc == 0 && TbVR == 0) {
	crosspoint = crossing(DeltaAltoud, DeltaAlt, MinTAVoud, MinTAVact);
	Tc = TimePointer - TimeStep * (1 - crosspoint);
      }
      if (DeltaAlt < MinTAVact && Ta == 0 && Tc != 0) {
	crosspoint = crossing(DeltaAltoud, DeltaAlt, MinTAVoud, MinTAVact);
	Ta = TimePointer - TimeStep * (1 - crosspoint);
      }
    } while (fabs(TimePointer - RiseSetS) <= MaxTryHours / 24.0 && Ta == 0 && !((TbVR != 0 && (TypeEvent == 3 || TypeEvent == 4) && (strncmp(ObjectName, "moon", 4) != 0 && strncmp(ObjectName, "venus", 5) != 0 && strncmp(ObjectName, "mercury", 7) != 0))));
    if (RS == 2) {
      TfirstVR = Tc;
      TlastVR = Ta;
    } else {
      TfirstVR = Ta;
      TlastVR = Tc;
    }
    if (TfirstVR == 0 && TlastVR == 0) {
      if (RS == 1)
        TfirstVR = TbVR - 0.000001;
      else
        TlastVR = TbVR + 0.000001;
    }
    if (!noriseO) {
      if (RS == 1)
        TfirstVR = mymax(TfirstVR, RiseSetO);
      else
        TlastVR = mymin(TlastVR, RiseSetO);
    }
    TvisVR = TJD_INVALID;
    if (TlastVR != 0 && TfirstVR != 0)
      TvisVR = TlastVR - TfirstVR;
    if (TlastVR == 0) TlastVR = TJD_INVALID;
    if (TbVR == 0) TbVR = TJD_INVALID;
    if (TfirstVR == 0) TfirstVR = TJD_INVALID;
output_heliacal_pheno:

  darr[0] = AltO;
  darr[1] = AppAltO;
  darr[2] = GeoAltO;
  darr[3] = AziO;
  darr[4] = AltS;
  darr[5] = AziS;
  darr[6] = TAVact;
  darr[7] = ARCVact;
  darr[8] = DAZact;
  darr[9] = ARCLact;
  darr[10] = kact;
  darr[11] = MinTAV;
  darr[12] = TfirstVR;
  darr[13] = TbVR;
  darr[14] = TlastVR;
  darr[15] = TbYallop;
  darr[16] = WMoon;
  darr[17] = qYal;
  darr[18] = qCrit;
  darr[19] = ParO;
  darr[20] = MagnO;
  darr[21] = RiseSetO;
  darr[22] = RiseSetS;
  darr[23] = Lag;
  darr[24] = TvisVR;
  darr[25] = LMoon;
  darr[26] = elong;
  darr[27] = illum;
  return OK;
}

#if 0
int32 HeliacalJDut(double JDNDaysUTStart, double Age, double SN, double Lat, double Longitude, double HeightEye, double Temperature, double Pressure, double RH, double VR, char *ObjectName, int TypeEvent, char *AVkind, double *dret, char *serr)
{
  double dgeo[3], datm[4], dobs[6];
  int32 helflag = SE_HELFLAG_HIGH_PRECISION;
  helflag |= SE_HELFLAG_AVKIND_VR;
  dgeo[0] = Longitude;
  dgeo[1] = Lat;
  dgeo[2] = HeightEye;
  datm[0] = Pressure;
  datm[1] = Temperature;
  datm[2] = RH;
  datm[3] = VR;
  dobs[0] = Age;
  dobs[1] = SN;
  return swe_heliacal_ut(JDNDaysUTStart, dgeo, datm, dobs, ObjectName, TypeEvent, helflag, dret, serr);
}
#endif

static double get_synodic_period(int Planet)
{

  switch(Planet) {
    case SE_MOON: return 29.530588853;
    case SE_MERCURY: return 115.8775;
    case SE_VENUS: return 583.9214;
    case SE_MARS: return 779.9361;
    case SE_JUPITER: return 398.8840;
    case SE_SATURN: return 378.0919;
    case SE_URANUS: return 369.6560;
    case SE_NEPTUNE: return 367.4867;
    case SE_PLUTO: return 366.7207;
  }
  return 366;
}

static int32 moon_event_arc_vis(double JDNDaysUTStart, double *dgeo, double *datm, double *dobs, int32 TypeEvent, int32 helflag, double *dret, char *serr)
{
  double x[20], MinTAV, MinTAVoud, OldestMinTAV;
  double phase1, phase2, JDNDaysUT, JDNDaysUTi;
  double tjd_moonevent, tjd_moonevent_start;
  double DeltaAltoud, TimeCheck, LocalminCheck;
  double AltS, AltO, DeltaAlt = 90;
  char ObjectName[30];
  int32 iflag, Daystep, goingup, Planet, retval;
  int32 avkind = helflag & SE_HELFLAG_AVKIND;
  int32 epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  dret[0] = JDNDaysUTStart;
  if (avkind == 0)
    avkind = SE_HELFLAG_AVKIND_VR;
  if (avkind != SE_HELFLAG_AVKIND_VR) {
    if (serr != NULL)
      strcpy(serr, "error: in valid AV kind for the moon");
    return ERR;
  }
  if (TypeEvent == 1 || TypeEvent == 2) {
    if (serr != NULL)
      strcpy(serr, "error: the moon has no morning first or evening last");
    return ERR;
  }
  strcpy(ObjectName, "moon");
  Planet = SE_MOON;
  iflag = SEFLG_TOPOCTR | SEFLG_EQUATORIAL | epheflag;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT|SEFLG_TRUEPOS;
  Daystep = 1;
  if (TypeEvent == 3) {

    TypeEvent = 2;
  } else {

    TypeEvent = 1;
    Daystep = -Daystep;
  }

  JDNDaysUT = JDNDaysUTStart;

  if (TypeEvent == 1) JDNDaysUT = JDNDaysUT + 30;

  swe_pheno_ut(JDNDaysUT, Planet, iflag, x, serr);
  phase2 = x[0];
  goingup = 0;
  do {
    JDNDaysUT = JDNDaysUT + Daystep;
    phase1 = phase2;
    swe_pheno_ut(JDNDaysUT, Planet, iflag, x, serr);
    phase2 = x[0];
    if (phase2 > phase1)
      goingup = 1;
  } while (goingup == 0 || (goingup == 1 && (phase2 > phase1)));

  JDNDaysUT = JDNDaysUT - Daystep;

  JDNDaysUTi = JDNDaysUT;
  JDNDaysUT = JDNDaysUT - Daystep;
  MinTAVoud = 199;
  do {
    JDNDaysUT = JDNDaysUT + Daystep;
    if ((retval = RiseSet(JDNDaysUT, dgeo, datm, ObjectName, TypeEvent, helflag, 0, &tjd_moonevent, serr)) != OK)
      return retval;
    tjd_moonevent_start = tjd_moonevent;
    MinTAV = 199;
    OldestMinTAV = MinTAV;
    do {
      OldestMinTAV = MinTAVoud;
      MinTAVoud = MinTAV;
      DeltaAltoud = DeltaAlt;
      tjd_moonevent = tjd_moonevent - 1.0 / 60.0 / 24.0 * Sgn(Daystep);
      if (ObjectLoc(tjd_moonevent, dgeo, datm, "sun", 0, helflag, &AltS, serr) == ERR)
	return ERR;
      if (ObjectLoc(tjd_moonevent, dgeo, datm, ObjectName, 0, helflag, &AltO, serr) == ERR)
	return ERR;
      DeltaAlt = AltO - AltS;
      if (DeterTAV(dobs, tjd_moonevent, dgeo, datm, ObjectName, helflag, &MinTAV, serr) == ERR)
        return ERR;
      TimeCheck = tjd_moonevent - LocalMinStep / 60.0 / 24.0 * Sgn(Daystep);
      if (DeterTAV(dobs, TimeCheck, dgeo, datm, ObjectName, helflag, &LocalminCheck, serr) == ERR)
        return ERR;

    } while ((MinTAV <= MinTAVoud || LocalminCheck < MinTAV) && fabs(tjd_moonevent - tjd_moonevent_start) < 120.0 / 60.0 / 24.0);

  } while (DeltaAltoud < MinTAVoud && fabs(JDNDaysUT - JDNDaysUTi) < 15);
  if (fabs(JDNDaysUT - JDNDaysUTi) < 15) {
    tjd_moonevent += (1 - x2min(MinTAV, MinTAVoud, OldestMinTAV)) * Sgn(Daystep) / 60.0 / 24.0;
  } else {
    strcpy(serr, "no date found for lunar event");
    return ERR;
  }
  dret[0] = tjd_moonevent;
  return OK;
}

static int32 heliacal_ut_arc_vis(double JDNDaysUTStart, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 TypeEventIn, int32 helflag, double *dret, char *serr_ret)
{
  double x[6];
  double xin[2];
  double xaz[2];
  double dang[3];
  double objectmagn = 0, maxlength, DayStep;
  double JDNDaysUT, JDNDaysUTfinal, JDNDaysUTstep, JDNDaysUTstepoud, JDNarcvisUT, tjd_tt, tret, OudeDatum, JDNDaysUTinp = JDNDaysUTStart, JDNDaysUTtijd;
  double ArcusVis, ArcusVisDelta, ArcusVisPto, ArcusVisDeltaoud;
  double Trise, sunsangle, Theliacal, Tdelta, Angle;
  double TimeStep, TimePointer, OldestMinTAV, MinTAVoud, MinTAVact, extrax, TbVR = 0;
  double AziS, AltS, AziO, AltO, DeltaAlt;
  double direct, Pressure, Temperature, d;
  int32 epheflag, retval = OK;
  int32 iflag, Planet, eventtype;
  int32 TypeEvent = TypeEventIn;
  int doneoneday;
  char serr[AS_MAXCH];
  *dret = JDNDaysUTStart;
  *serr = '\0';
  Planet = DeterObject(ObjectName);
  Pressure = datm[0];
  Temperature = datm[1];

  if ((retval = Magnitude(JDNDaysUTStart, dgeo, ObjectName, helflag, &objectmagn, serr)) == ERR)
    goto swe_heliacal_err;
  epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  iflag = SEFLG_TOPOCTR | SEFLG_EQUATORIAL | epheflag;
  if (!(helflag & SE_HELFLAG_HIGH_PRECISION))
    iflag |= SEFLG_NONUT | SEFLG_TRUEPOS;

  switch(Planet) {
    case SE_MERCURY:
      DayStep = 1; maxlength = 100; break;
    case SE_VENUS:
      DayStep = 64; maxlength = 384; break;
    case SE_MARS:
      DayStep = 128; maxlength = 640; break;
    case SE_JUPITER:
      DayStep = 64; maxlength = 384; break;
    case SE_SATURN:
      DayStep = 64; maxlength = 256; break;
    default:
      DayStep = 64; maxlength = 256; break;
  }

  eventtype = TypeEvent;
  if (eventtype == 2) DayStep = -DayStep;

  if (eventtype == 4) {
      eventtype = 1;
      DayStep = -DayStep;
  }

  if (eventtype == 3) eventtype = 2;
  eventtype |= SE_BIT_DISC_CENTER;

  {

    JDNDaysUT = JDNDaysUTStart;

    JDNDaysUTfinal = JDNDaysUT + maxlength;
    JDNDaysUT = JDNDaysUT - 1;
    if (DayStep < 0) {
      JDNDaysUTtijd = JDNDaysUT;
      JDNDaysUT = JDNDaysUTfinal;
      JDNDaysUTfinal = JDNDaysUTtijd;
    }

    JDNDaysUTstep = JDNDaysUT - DayStep;
    doneoneday = 0;
    ArcusVisDelta = 199;
    ArcusVisPto = -5.55;
    do {
      if (fabs(DayStep) == 1) doneoneday = 1;
      do {

	JDNDaysUTstepoud = JDNDaysUTstep;
	ArcusVisDeltaoud = ArcusVisDelta;
	JDNDaysUTstep = JDNDaysUTstep + DayStep;

	if ((retval = my_rise_trans(JDNDaysUTstep, SE_SUN, "", eventtype, helflag, dgeo, datm, &tret, serr)) == ERR)
	  goto swe_heliacal_err;

	tjd_tt = tret + swe_deltat_ex(tret, epheflag, serr);
	if ((retval = swe_calc(tjd_tt, SE_SUN, iflag, x, serr)) == ERR)
	  goto swe_heliacal_err;
	xin[0] = x[0];
	xin[1] = x[1];
	swe_azalt(tret, SE_EQU2HOR, dgeo, Pressure, Temperature, xin, xaz);
	Trise = HourAngle(xaz[1], x[1], dgeo[1]);
	sunsangle = ArcusVisPto;
	if (helflag & SE_HELFLAG_AVKIND_MIN7) sunsangle = -7;
	if (helflag & SE_HELFLAG_AVKIND_MIN9) sunsangle = -9;
	Theliacal = HourAngle(sunsangle, x[1], dgeo[1]);
	Tdelta = Theliacal - Trise;
	if (TypeEvent == 2 || TypeEvent== 3) Tdelta = -Tdelta;

	JDNarcvisUT = tret - Tdelta / 24;
	tjd_tt = JDNarcvisUT + swe_deltat_ex(JDNarcvisUT, epheflag, serr);

	if ((retval = swe_calc(tjd_tt, SE_SUN, iflag, x, serr)) == ERR)
	  goto swe_heliacal_err;
	xin[0] = x[0];
	xin[1] = x[1];
	swe_azalt(JDNarcvisUT, SE_EQU2HOR, dgeo, Pressure, Temperature, xin, xaz);
	AziS = xaz[0] + 180;
	if (AziS >= 360) AziS = AziS - 360;
	AltS = xaz[1];

#if 0
  double AltM, AziM;
	if ((retval = swe_calc(tjd_tt, SE_MOON, iflag, x, serr)) == ERR)
	  goto swe_heliacal_err;
	xin[0] = x[0];
	xin[1] = x[1];
	swe_azalt(JDNarcvisUT, SE_EQU2HOR, dgeo, Pressure, Temperature, xin, xaz);
	AziM = xaz[0] + 180;
	if (AziM >= 360) AziM = AziM - 360;
	AltM = xaz[1];
#endif

	if (Planet != -1) {
	  if ((retval = swe_calc(tjd_tt, Planet, iflag, x, serr)) == ERR)
	    goto swe_heliacal_err;

	  if ((retval = Magnitude(JDNarcvisUT, dgeo, ObjectName, helflag, &objectmagn, serr)) == ERR)
	    goto swe_heliacal_err;
	} else {
	  if ((retval = call_swe_fixstar(ObjectName, tjd_tt, iflag, x, serr)) == ERR)
	    goto swe_heliacal_err;
	}
	xin[0] = x[0];
	xin[1] = x[1];
	swe_azalt(JDNarcvisUT, SE_EQU2HOR, dgeo, Pressure, Temperature, xin, xaz);
	AziO = xaz[0] + 180;
	if (AziO >= 360) AziO = AziO - 360;
	AltO = xaz[1];

	DeltaAlt = AltO - AltS;

	if ((retval = HeliacalAngle(objectmagn, dobs, AziO, -1, 0, JDNarcvisUT, AziS, dgeo, datm, helflag, dang, serr)) == ERR)
	  goto swe_heliacal_err;
	ArcusVis = dang[1];
	ArcusVisPto = dang[2];
	ArcusVisDelta = DeltaAlt - ArcusVis;

      } while ((ArcusVisDeltaoud > 0 || ArcusVisDelta < 0) && (JDNDaysUTfinal - JDNDaysUTstep) * Sgn(DayStep) > 0);
      if (doneoneday == 0 && (JDNDaysUTfinal - JDNDaysUTstep) * Sgn(DayStep) > 0) {

	ArcusVisDelta = ArcusVisDeltaoud;
	DayStep = ((int) (fabs(DayStep) / 2.0)) * Sgn(DayStep);
	JDNDaysUTstep = JDNDaysUTstepoud;
      }
    } while (doneoneday == 0 && (JDNDaysUTfinal - JDNDaysUTstep) * Sgn(DayStep) > 0);
  }
  d = (JDNDaysUTfinal - JDNDaysUTstep) * Sgn(DayStep);
  if (d <= 0 || d >= maxlength) {
    dret[0] = JDNDaysUTinp;
    retval = -2;
    sprintf(serr, "heliacal event not found within maxlength %f\n", maxlength);
    goto swe_heliacal_err;
  }
#if 0
  if (helflag & SE_HELFLAG_AVKIND_VR) {
    double darr[40];
    if (swe_heliacal_pheno_ut(JDNarcvisUT, dgeo, datm, dobs, ObjectName, TypeEvent, helflag, darr, serr) != OK)
      return ERR;
    JDNarcvisUT = darr[13];
    }
  }
#endif
  direct = TimeStepDefault / 24.0 / 60.0;
  if (DayStep < 0) direct = -direct;
  if (helflag & SE_HELFLAG_AVKIND_VR) {

    TimeStep = direct;
    TbVR = 0;
    TimePointer = JDNarcvisUT;
    if (DeterTAV(dobs, TimePointer, dgeo, datm, ObjectName, helflag, &OldestMinTAV, serr) == ERR)
      return ERR;
    TimePointer = TimePointer + TimeStep;
    if (DeterTAV(dobs, TimePointer, dgeo, datm, ObjectName, helflag, &MinTAVoud, serr) == ERR)
      return ERR;
    if (MinTAVoud > OldestMinTAV) {
      TimePointer = JDNarcvisUT;
      TimeStep = -TimeStep;
      MinTAVact = OldestMinTAV;
    } else {
      MinTAVact = MinTAVoud;
      MinTAVoud = OldestMinTAV;
    }

    do {
      TimePointer = TimePointer + TimeStep;
      OldestMinTAV = MinTAVoud;
      MinTAVoud = MinTAVact;
      if (DeterTAV(dobs, TimePointer, dgeo, datm, ObjectName, helflag, &MinTAVact, serr) == ERR)
        return ERR;
      if (MinTAVoud < MinTAVact) {
	extrax = x2min(MinTAVact, MinTAVoud, OldestMinTAV);
	TbVR = TimePointer - (1 - extrax) * TimeStep;
      }
    } while (TbVR == 0);
    JDNarcvisUT = TbVR;
  }

  if (helflag & SE_HELFLAG_AVKIND_PTO) {
    do {
      OudeDatum = JDNarcvisUT;
      JDNarcvisUT = JDNarcvisUT - direct;
      tjd_tt = JDNarcvisUT + swe_deltat_ex(JDNarcvisUT, epheflag, serr);
      if (Planet != -1) {
	if ((retval = swe_calc(tjd_tt, Planet, iflag, x, serr)) == ERR)
	  goto swe_heliacal_err;
      } else {
	if ((retval = call_swe_fixstar(ObjectName, tjd_tt, iflag, x, serr)) == ERR)
	  goto swe_heliacal_err;
      }
      xin[0] = x[0];
      xin[1] = x[1];
      swe_azalt(JDNarcvisUT, SE_EQU2HOR, dgeo, Pressure, Temperature, xin, xaz);
      Angle = xaz[1];
    } while (Angle > 0);
    JDNarcvisUT = (JDNarcvisUT + OudeDatum) / 2.0;
  }
  if (JDNarcvisUT < -9999999 || JDNarcvisUT > 9999999) {
    dret[0] = JDNDaysUT;
    strcpy(serr, "no heliacal date found");
    retval = ERR;
    goto swe_heliacal_err;
  }
  dret[0] = JDNarcvisUT;
swe_heliacal_err:
  if (serr_ret != NULL && *serr != '\0')
    strcpy(serr_ret, serr);
  return retval;
}

static int32 get_asc_obl(double tjd, int32 ipl, char *star, int32 iflag, double
 *dgeo, AS_BOOL desc_obl, double *daop, char *serr)
{
  int32 retval;
  int32 epheflag = iflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  double x[6], adp;
  char s[AS_MAXCH];
  char star2[AS_MAXCH];
  strcpy(star2, star);
  if (ipl == -1) {
    if ((retval = swe_fixstar(star2, tjd, epheflag | SEFLG_EQUATORIAL, x, serr)) == ERR)
      return ERR;
  } else {
    if ((retval = swe_calc(tjd, ipl, epheflag | SEFLG_EQUATORIAL, x, serr)) == ERR)
      return ERR;
  }
  adp = tan(dgeo[1] * DEGTORAD) * tan(x[1] * DEGTORAD);
  if (fabs(adp) > 1) {
    if (star != NULL && *star != '\0')
      strcpy(s, star);
    else
      swe_get_planet_name(ipl, s);
    sprintf(serr, "%s is circumpolar, cannot calculate heliacal event", s);
    return -2;
  }
  adp = asin(adp) / DEGTORAD;
  if (desc_obl)
    *daop = x[0] + adp;
  else
    *daop = x[0] - adp;
  *daop = swe_degnorm(*daop);
  return OK;
}

#if 0
static int32 get_asc_obl_old(double tjd, int32 ipl, char *star, int32 iflag, double *dgeo, AS_BOOL desc_obl, double *daop, char *serr)
{
  int32 retval;
  int32 epheflag = iflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  double x[6], adp;
  char s[AS_MAXCH];
  if (star != NULL && *star != '\0') {
    if ((retval = call_swe_fixstar(star, tjd, epheflag | SEFLG_EQUATORIAL, x, serr)) == ERR)
      return ERR;
  } else {
    if ((retval = swe_calc(tjd, ipl, epheflag | SEFLG_EQUATORIAL, x, serr)) == ERR)
      return ERR;
  }
  adp = tan(dgeo[1] * DEGTORAD) * tan(x[1] * DEGTORAD);
  if (fabs(adp) > 1) {
    if (star != NULL && *star != '\0')
      strcpy(s, star);
    else
      swe_get_planet_name(ipl, s);
    sprintf(serr, "%s is circumpolar, cannot calculate heliacal event", s);
    return -2;
  }
  adp = asin(adp) / DEGTORAD;
  if (desc_obl)
    *daop = x[0] + adp;
  else
    *daop = x[0] - adp;
  *daop = swe_degnorm(*daop);
  return OK;
}
#endif

static int32 get_asc_obl_diff(double tjd, int32 ipl, char *star, int32 iflag, double *dgeo, AS_BOOL desc_obl, AS_BOOL is_acronychal, double *dsunpl, char *serr)
{
  int32 retval = OK;
  double aosun, aopl;

  retval = get_asc_obl(tjd, SE_SUN, "", iflag, dgeo, desc_obl, &aosun, serr);
  if (retval != OK)
    return retval;
  if (is_acronychal) {
    if (desc_obl == TRUE)
      desc_obl = FALSE;
    else
      desc_obl = TRUE;
  }

  retval = get_asc_obl(tjd, ipl, star, iflag, dgeo, desc_obl, &aopl, serr);
  if (retval != OK)
    return retval;
  *dsunpl = swe_degnorm(aosun - aopl);
  if (is_acronychal)
    *dsunpl = swe_degnorm(*dsunpl - 180);
  if (*dsunpl > 180) *dsunpl -= 360;
  return OK;
}

#if 0
static int32 get_asc_obl_diff_old(double tjd, int32 ipl, char *star, int32 iflag, double *dgeo, AS_BOOL desc_obl, double *dsunpl, char *serr)
{
  int32 retval = OK;
  double aosun, aopl;

  retval = get_asc_obl(tjd, SE_SUN, "", iflag, dgeo, desc_obl, &aosun, serr);
  if (retval != OK)
    return retval;

  retval = get_asc_obl(tjd, ipl, star, iflag, dgeo, desc_obl, &aopl, serr);
  if (retval != OK)
    return retval;
  *dsunpl = swe_degnorm(aosun - aopl);
  return OK;
}
#endif

static const double tcon[] =
{
  0, 0,
  2451550, 2451550,
  2451604, 2451670,
  2451980, 2452280,
  2451727, 2452074,
  2451673, 2451877,
  2451675, 2451868,
  2451581, 2451768,
  2451568, 2451753,
};

static int32 find_conjunct_sun(double tjd_start, int32 ipl, int32 helflag, int32 TypeEvent, double *tjd, char *serr)
{
  int32 epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  int i;
  double tjdcon, tjd0, ds, dsynperiod, x[6], xs[6], daspect = 0;
  if (ipl >= SE_MARS && TypeEvent >= 3)
    daspect = 180;
  i = (TypeEvent - 1) / 2 + ipl * 2;
  tjd0 = tcon[i];
  dsynperiod = get_synodic_period(ipl);
  tjdcon = tjd0 + ((floor) ((tjd_start - tjd0) / dsynperiod) + 1) * dsynperiod;
  ds = 100;
  while (ds > 0.5) {
    if (swe_calc(tjdcon, ipl, epheflag|SEFLG_SPEED, x, serr) == ERR)
      return ERR;
    if (swe_calc(tjdcon, SE_SUN, epheflag|SEFLG_SPEED, xs, serr) == ERR)
      return ERR;
    ds = swe_degnorm(x[0] - xs[0] - daspect);
    if (ds > 180) ds -= 360;
    tjdcon -= ds / (x[3] - xs[3]);
  }
  *tjd = tjdcon;
  return OK;
}

static int32 get_asc_obl_with_sun(double tjd_start, int32 ipl, char *star, int32 helflag, int32 evtyp, double dperiod, double *dgeo, double *tjdret, char *serr)
{
  int i, retval;
  int32 is_acronychal = FALSE;
  int32 epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  double dsunpl = 1, dsunpl_save, dsunpl_test, tjd, daystep;
  AS_BOOL desc_obl = FALSE, retro = FALSE;
  if (evtyp == SE_EVENING_LAST || evtyp == SE_EVENING_FIRST)
    desc_obl = TRUE;
  if (evtyp == SE_MORNING_FIRST || evtyp == SE_EVENING_LAST)
    retro = TRUE;
  if (evtyp == SE_ACRONYCHAL_RISING)
    desc_obl = TRUE;
  if (evtyp == SE_ACRONYCHAL_RISING || evtyp ==  SE_ACRONYCHAL_SETTING) {
    is_acronychal = TRUE;
    if (ipl != SE_MOON)
      retro = TRUE;
  }

  tjd = tjd_start;
  dsunpl_save = -999999999;
  retval = get_asc_obl_diff(tjd, ipl, star, epheflag, dgeo, desc_obl, is_acronychal, &dsunpl, serr);
  if (retval != OK)
    return retval;
  daystep = 20;
  i = 0;
  while (dsunpl_save == -999999999 ||

      fabs(dsunpl) + fabs(dsunpl_save) > 180 ||
      (retro && !(dsunpl_save < 0 && dsunpl >= 0)) ||
      (!retro && !(dsunpl_save >= 0 && dsunpl < 0))) {
    i++;
    if (i > 5000) {
      sprintf(serr, "loop in get_asc_obl_with_sun() (1)");
      return ERR;
    }
    dsunpl_save = dsunpl;
    tjd += 10.0;
    if (dperiod > 0 && tjd - tjd_start > dperiod)
      return -2;
    retval = get_asc_obl_diff(tjd, ipl, star, epheflag, dgeo, desc_obl, is_acronychal, &dsunpl, serr);
    if (retval != OK)
      return retval;
  }
  tjd_start = tjd - daystep;
  daystep /= 2.0;
  tjd = tjd_start + daystep;
  retval = get_asc_obl_diff(tjd, ipl, star, epheflag, dgeo, desc_obl, is_acronychal, &dsunpl_test, serr);
  if (retval != OK)
    return retval;
  i = 0;
  while (fabs(dsunpl) > 0.00001) {
    i++;
    if (i > 5000) {
      sprintf(serr, "loop in get_asc_obl_with_sun() (2)");
      return ERR;
    }
    if (dsunpl_save * dsunpl_test >= 0) {
      dsunpl_save = dsunpl_test;
      tjd_start = tjd;
    } else {
      dsunpl = dsunpl_test;
    }
    daystep /= 2.0;
    tjd = tjd_start + daystep;
    retval = get_asc_obl_diff(tjd, ipl, star, epheflag, dgeo, desc_obl, is_acronychal, &dsunpl_test, serr);
    if (retval != OK)
      return retval;
  }
  *tjdret = tjd;
  return OK;
}

#if 0

static int32 get_asc_obl_with_sun_old(double tjd_start, int32 ipl, char *star, int32 helflag, int32 TypeEvent, double *dgeo, double *tjdret, char *serr)
{
  int retval;
  int32 epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  double dsunpl = 1, tjd, daystep, dsunpl_save;
  double dsynperiod = 367;
  double dangsearch = 0;
  AS_BOOL desc_obl = FALSE;
  if (TypeEvent == 2 || TypeEvent == 3)
    desc_obl = TRUE;
  if (TypeEvent == 3 || TypeEvent == 4)
    dangsearch = 180;

  daystep = dsynperiod;
  tjd = tjd_start;
  retval = get_asc_obl_diff(tjd, ipl, star, epheflag, dgeo, desc_obl, &dsunpl, serr);
  if (retval != OK)
    return retval;
  while (dsunpl < 359.99999) {
    dsunpl_save = dsunpl;
    daystep /= 2.0;
    retval = get_asc_obl_diff(tjd + daystep, ipl, star, epheflag, dgeo, desc_obl, &dsunpl, serr);
    if (retval != OK)
      return retval;
    if (dsunpl > dsunpl_save)
      tjd += daystep;
    else
      dsunpl = dsunpl_save;
  }
  *tjdret = tjd;
  return OK;
}
#endif

#if 0

static int32 get_asc_obl_acronychal(double tjd_start, int32 ipl, char *star, int32 helflag, int32 TypeEvent, double *dgeo, double *tjdret, char *serr)
{
  int retval;
  int32 epheflag = helflag & (SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH);
  double dsunpl = 1, tjd, daystep, dsunpl_save;
  double dsynperiod = 367;
  double aosun, aopl;
  AS_BOOL sun_desc = TRUE, obj_desc = FALSE;
  daystep = dsynperiod;
  tjd = tjd_start;
  if (TypeEvent == 4) {
    sun_desc = FALSE;
    obj_desc = TRUE;
  }

  retval = get_asc_obl(tjd, SE_SUN, "", epheflag, dgeo, sun_desc, &aosun, serr);
  if (retval != OK)
    return retval;

  retval = get_asc_obl(tjd, ipl, star, epheflag, dgeo, obj_desc, &aopl, serr);
  if (retval != OK)
    return retval;
  dsunpl = swe_degnorm(aosun - aopl + 180);
  while (dsunpl < 359.99999) {
    dsunpl_save = dsunpl;
    daystep /= 2.0;

    retval = get_asc_obl(tjd+daystep, SE_SUN, "", epheflag, dgeo, sun_desc, &aosun, serr);
    if (retval != OK)
      return retval;

    retval = get_asc_obl(tjd+daystep, ipl, star, epheflag, dgeo, obj_desc, &aopl, serr);
    if (retval != OK)
      return retval;
    dsunpl = swe_degnorm(aosun - aopl + 180);
    if (dsunpl > dsunpl_save)
      tjd += daystep;
    else
      dsunpl = dsunpl_save;
  }
  *tjdret = tjd;
  return OK;
}
#endif

static int32 get_heliacal_day(double tjd, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 helflag, int32 TypeEvent, double *thel, char *serr)
{
  int32 i, visible_at_sunsetrise, is_rise_or_set = 0, ndays, retval, retval_old;
  double direct_day = 0, direct_time = 0, tfac, tend, daystep, tday, vdelta, tret;
  double darr[30], vd, dmag, div;
  int32 ipl = DeterObject(ObjectName);

  switch (TypeEvent) {

    case 1: is_rise_or_set = SE_CALC_RISE;
      direct_day = 1; direct_time = -1;
      break;

    case 2: is_rise_or_set = SE_CALC_SET;
      direct_day = -1; direct_time = 1;
      break;

    case 3: is_rise_or_set = SE_CALC_SET;
      direct_day = 1; direct_time = 1;
      break;

    case 4: is_rise_or_set = SE_CALC_RISE;
      direct_day = -1; direct_time = -1;
      break;
  }
  tfac = 1;
  switch (ipl) {
    case SE_MOON:
      ndays = 16;
      daystep = 1;
      break;
    case SE_MERCURY:
      ndays = 60; tjd -= 0 * direct_day;
      daystep = 5;
      tfac = 5;
      break;
    case SE_VENUS:
      ndays = 300; tjd -= 30 * direct_day;
      daystep = 5;
      if (TypeEvent >= 3) {
	daystep = 15;
        tfac = 3;
      }
      break;
    case SE_MARS:
      ndays = 400;
      daystep = 15;
      tfac = 5;
      break;
    case SE_SATURN:
      ndays = 300;
      daystep = 20;
      tfac = 5;
      break;
    case -1:
      ndays = 300;
      if (call_swe_fixstar_mag(ObjectName, &dmag, serr) == ERR) {
	return ERR;
      }
      daystep = 15;
      tfac = 10;
      if (dmag > 2) {
        daystep = 15;
      }
      if (dmag < 0) {
	tfac = 3;
      }
      break;
    default:
      ndays = 300;
      daystep = 15;
      tfac = 3;
      break;
  }
  tend = tjd + ndays * direct_day;
  retval_old = -2;
  for (tday = tjd, i = 0;
       (direct_day > 0 && tday < tend) || (direct_day < 0 && tday > tend);
       tday += daystep * direct_day, i++) {
    vdelta = -100;
    if (i > 0)
      tday -= 0.3 * direct_day;
    if ((retval = my_rise_trans(tday, SE_SUN, "", is_rise_or_set, helflag, dgeo, datm, &tret, serr)) == ERR) {
      return ERR;
    }

    if (retval == -2) {
      retval_old = retval;
      continue;
    }
    retval = swe_vis_limit_mag(tret, dgeo, datm, dobs, ObjectName, helflag, darr, serr);
    if (retval == ERR)
      return ERR;
#if 1

    if (retval_old == -2 && retval >= 0 && daystep > 1) {
      retval_old = retval;
      tday -= daystep * direct_day;
      daystep = 1;

      if (ipl >= SE_MARS || ipl == -1)
        daystep = 5;
      continue;
    }
    retval_old = retval;
#endif

    if (retval == -2)
      continue;
    vdelta = darr[0] - darr[7];

    div = 1440.0;

    vd = -1;
    visible_at_sunsetrise = 1;
    while (retval != -2 && (vd = darr[0] - darr[7]) < 0) {
      visible_at_sunsetrise = 0;
      if (vd < -1.0)
	tret += 5.0 / div * direct_time * tfac;
      else if (vd < -0.5)
	tret += 2.0 / div * direct_time * tfac;
      else if (vd < -0.1)
	tret += 1.0 / div * direct_time * tfac;
      else
	tret += 1.0 / div * direct_time;
      retval = swe_vis_limit_mag(tret, dgeo, datm, dobs, ObjectName, helflag, darr, serr);
      if (retval == ERR)
	return ERR;
    }

    if (visible_at_sunsetrise) {
      for (i = 0; i < 10; i++) {
	if ((retval = swe_vis_limit_mag(tret + 1.0 / div * direct_time, dgeo, datm, dobs, ObjectName, helflag, darr, serr)) >= 0
	  && darr[0] - darr[7] > vd) {
          vd = darr[0] - darr[7];
	  tret += 1.0 / div * direct_time;
	}
      }
    }
    vdelta = darr[0] - darr[7];

    if (vdelta > 0) {
      if ((ipl >= SE_MARS || ipl == -1) && daystep > 1) {
	tday -= daystep * direct_day;
	daystep = 1;
      } else {
	*thel = tret;
	return OK;
      }
    }
  }
  sprintf(serr, "heliacal event does not happen");

  return -2;
}

static int32 time_optimum_visibility(double tjd, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 helflag, double *tret, char *serr)
{
  int32 retval, retval_sv, i;
  double t1, t2, vl1, vl2, d, darr[10], phot_scot_opic, phot_scot_opic_sv;

  int t_has_changed;
  *tret = tjd;
  retval = swe_vis_limit_mag(tjd, dgeo, datm, dobs, ObjectName, helflag, darr, serr);
  if (retval == ERR) return ERR;
  retval_sv = retval;

  t1 = tjd;
  t2 = tjd;
  vl1 = -1;
  vl2 = -1;

  phot_scot_opic_sv = retval & SE_SCOTOPIC_FLAG;
  for (i = 0, d = 100.0 / 86400.0; i < 3; i++, d /= 10.0) {

    t1 += d;
    t_has_changed = 0;
    while((retval = swe_vis_limit_mag(t1 - d, dgeo, datm, dobs, ObjectName, helflag, darr, serr)) >= 0
	&& darr[0] > darr[7]
	&& darr[0] - darr[7] > vl1) {
      t1 -= d; vl1 = darr[0] - darr[7];
      t_has_changed = 1;

      retval_sv = retval;
      phot_scot_opic_sv = retval & SE_SCOTOPIC_FLAG;

    }
    if (t_has_changed == 0)
      t1 -= d;
    if (retval == ERR) return ERR;
  }
  for (i = 0, d = 100.0 / 86400.0; i < 3; i++, d /= 10.0) {
    t2 -= d;
    t_has_changed = 0;
    while((retval = swe_vis_limit_mag(t2 + d, dgeo, datm, dobs, ObjectName, helflag, darr, serr)) >= 0
        && darr[0] > darr[7]
	&& darr[0] - darr[7] > vl2) {
      t2 += d; vl2 = darr[0] - darr[7];
      t_has_changed = 1;

      retval_sv = retval;
      phot_scot_opic_sv = retval & SE_SCOTOPIC_FLAG;

    }
    if (t_has_changed == 0)
      t2 += d;
    if (retval == ERR) return ERR;
  }
  if (vl2 > vl1)
    tjd = t2;
  else
    tjd = t1;

  *tret = tjd;

  if (retval >= 0) {

    phot_scot_opic = (retval & SE_SCOTOPIC_FLAG);
    if (phot_scot_opic_sv != phot_scot_opic) {

      printf ("hallo -2\n");
      return -2;
    }

    if (retval_sv & SE_MIXEDOPIC_FLAG) {
      printf ("hallo -2\n");
      return -2;
    }
  }
  return OK;
}

static int32 time_limit_invisible(double tjd, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 helflag, int32 direct, double *tret, char *serr)
{
  int32 retval, retval_sv, i, ncnt = 3;
  double d = 0, darr[10], phot_scot_opic, phot_scot_opic_sv;
  double d0 = 100.0 / 86400.0;
  *tret = tjd;
  if (strcmp(ObjectName, "moon") == 0) {
    d0 *= 10;
    ncnt = 4;
  }
  retval = swe_vis_limit_mag(tjd + d * direct, dgeo, datm, dobs, ObjectName, helflag, darr, serr);
  if (retval == ERR) return ERR;
  retval_sv = retval;
  phot_scot_opic_sv = retval & SE_SCOTOPIC_FLAG;
  for (i = 0, d = d0; i < ncnt; i++, d /= 10.0) {
    while((retval = swe_vis_limit_mag(tjd + d * direct, dgeo, datm, dobs, ObjectName, helflag, darr, serr)) >= 0
        && darr[0] > darr[7]) {
      tjd += d * direct;
      retval_sv = retval;
      phot_scot_opic_sv = retval & SE_SCOTOPIC_FLAG;

    }
  }

  *tret = tjd;

  *serr = '\0';
  if (retval >= 0) {

    phot_scot_opic = (retval & SE_SCOTOPIC_FLAG);
    if (phot_scot_opic_sv != phot_scot_opic) {

      return -2;
    }

    if (retval_sv & SE_MIXEDOPIC_FLAG) {
      return -2;
    }
  }
  return OK;
}

static int32 get_acronychal_day(double tjd, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 helflag, int32 TypeEvent, double *thel, char *serr) {
  double tret = 0, tret_dark = 0, darr[30], dtret;

  int32 retval, is_rise_or_set, direct;
  int32 ipl = DeterObject(ObjectName);
  helflag |= SE_HELFLAG_VISLIM_PHOTOPIC;

  if (TypeEvent == 3 || TypeEvent == 5) {
    is_rise_or_set = SE_CALC_RISE;

    direct = -1;
  } else {
    is_rise_or_set = SE_CALC_SET;

    direct = 1;
  }
  dtret = 999;
#if 0
  while (fabs(dtret) > 0.5) {
#else
  while (fabs(dtret) > 0.5 / 1440.0) {
#endif
    tjd += 0.7 * direct;
    if (direct < 0) tjd -= 1;
    retval = my_rise_trans(tjd, ipl, ObjectName, is_rise_or_set, helflag, dgeo, datm, &tjd, serr);
    if (retval == ERR) return ERR;
    retval = swe_vis_limit_mag(tjd, dgeo, datm, dobs, ObjectName, helflag, darr, serr);
    if (retval == ERR) return ERR;
    while(darr[0] < darr[7]) {
      tjd += 10.0 / 1440.0 * -direct;
      retval = swe_vis_limit_mag(tjd, dgeo, datm, dobs, ObjectName, helflag, darr, serr);
      if (retval == ERR) return ERR;
    }
    retval = time_limit_invisible(tjd, dgeo, datm, dobs, ObjectName, helflag | SE_HELFLAG_VISLIM_DARK, direct, &tret_dark, serr);
    if (retval == ERR) return ERR;
    retval = time_limit_invisible(tjd, dgeo, datm, dobs, ObjectName, helflag | SE_HELFLAG_VISLIM_NOMOON, direct, &tret, serr);
    if (retval == ERR) return ERR;
#if 0
    if (azalt_cart(tret_dark, dgeo, datm, ObjectName, helflag, darr, serr) == ERR)
      return ERR;
    if (azalt_cart(tret, dgeo, datm, ObjectName, helflag, darr+6, serr) == ERR)
      return ERR;
    dtret = acos(swi_dot_prod_unit(darr+3, darr+9)) / DEGTORAD;
#else
    dtret = fabs(tret - tret_dark);
#endif
  }
  if (azalt_cart(tret, dgeo, datm, "sun", helflag, darr, serr) == ERR)
    return ERR;
  *thel = tret;
  if (darr[1] < -12) {
    sprintf(serr, "acronychal rising/setting not available, %f", darr[1]);
    return OK;
  } else {
    sprintf(serr, "solar altitude, %f", darr[1]);
  }
  return OK;
}

static int32 get_heliacal_details(double tday, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 helflag, double *dret, char *serr)
{
  int32 i, retval, direct;
  AS_BOOL optimum_undefined, limit_1_undefined, limit_2_undefined;

  optimum_undefined = FALSE;
  retval = time_optimum_visibility(tday, dgeo, datm, dobs, ObjectName, helflag, &(dret[1]), serr);
  if (retval == ERR) return ERR;
  if (retval == -2) {
    retval = OK;
    optimum_undefined = TRUE;
  }

  direct = 1;
  if (TypeEvent == 1 || TypeEvent == 4)
    direct = -1;
  limit_1_undefined = FALSE;
  retval = time_limit_invisible(tday, dgeo, datm, dobs, ObjectName, helflag, direct, &(dret[0]), serr);
  if (retval == ERR) return ERR;
  if (retval == -2) {
    retval = OK;
    limit_1_undefined = TRUE;
  }

  direct *= -1;
  limit_2_undefined = FALSE;
  retval = time_limit_invisible(dret[1], dgeo, datm, dobs, ObjectName, helflag, direct, &(dret[2]), serr);
  if (retval == ERR) return ERR;
  if (retval == -2) {
    retval = OK;
    limit_2_undefined = TRUE;
  }

  if (TypeEvent == 2 || TypeEvent == 3) {
    tday = dret[2];
    dret[2] = dret[0];
    dret[0] = tday;
    i = (int) limit_1_undefined;
    limit_1_undefined = limit_2_undefined;
    limit_2_undefined = (AS_BOOL) i;
  }

  if (optimum_undefined || limit_1_undefined || limit_2_undefined) {
    sprintf(serr, "return values [");
    if (limit_1_undefined)
      strcat(serr, "0,");
    if (optimum_undefined)
      strcat(serr, "1,");
    if (limit_2_undefined)
      strcat(serr, "2,");
    strcat(serr, "] are uncertain due to change between photopic and scotopic vision");
  }
  return OK;
}

static int32 heliacal_ut_vis_lim(double tjd_start, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 TypeEventIn, int32 helflag, double *dret, char *serr_ret)
{
  int i;
  double d, darr[10], direct = 1, tjd, tday;
  int32 retval = OK, helflag2;
  int32 ipl;
  int32 TypeEvent = TypeEventIn;
  char serr[AS_MAXCH];
  for (i = 0; i < 10; i++)
    dret[i] = 0;
  *dret = tjd_start;
  *serr = '\0';
  ipl = DeterObject(ObjectName);
  if (ipl == SE_MERCURY)
    tjd = tjd_start - 30;
  else
    tjd = tjd_start - 50;

  helflag2 = helflag;

  if (ipl == SE_MERCURY || ipl == SE_VENUS || TypeEvent <= 2) {
    if (ipl == -1) {

    retval = get_asc_obl_with_sun(tjd, ipl, ObjectName, helflag, TypeEvent, 0, dgeo, &tjd, serr);
      if (retval != OK)
	goto swe_heliacal_err;
    } else {

      if ((retval = find_conjunct_sun(tjd, ipl, helflag, TypeEvent, &tjd, serr)) == ERR) {
	goto swe_heliacal_err;
      }
    }

    retval = get_heliacal_day(tjd, dgeo, datm, dobs, ObjectName, helflag2, TypeEvent, &tday, serr);
    if (retval != OK)
      goto swe_heliacal_err;

  } else {

    retval = get_asc_obl_with_sun(tjd, ipl, ObjectName, helflag, TypeEvent, 0, dgeo, &tjd, serr);
    if (retval != OK)
      goto swe_heliacal_err;
    tday = tjd;
    retval = get_acronychal_day(tjd, dgeo, datm, dobs, ObjectName, helflag2, TypeEvent, &tday, serr);
    if (retval != OK)
      goto swe_heliacal_err;
  }
  dret[0] = tday;
  if (!(helflag & SE_HELFLAG_NO_DETAILS)) {

    if (ipl == SE_MERCURY || ipl == SE_VENUS || TypeEvent <= 2) {
      retval = get_heliacal_details(tday, dgeo, datm, dobs, ObjectName, TypeEvent, helflag2, dret, serr);
      if (retval == ERR) goto swe_heliacal_err;
    } else if ((0)) {
      if (TypeEvent == 4 || TypeEvent == 6) direct = -1;
      for (i = 0, d = 100.0 / 86400.0; i < 3; i++, d /= 10.0) {
	while((retval = swe_vis_limit_mag(*dret + d * direct, dgeo, datm, dobs, ObjectName, helflag, darr, serr)) == -2 || (retval >= 0 && darr[0] < darr[7])) {
	  *dret += d * direct;
	}
      }

      if (retval == OK)
	*dret += 1.0 / 86400.0 * direct;
    }
  }
swe_heliacal_err:
  if (serr_ret != NULL && *serr != '\0')
    strcpy(serr_ret, serr);
  return retval;
}

static int32 moon_event_vis_lim(double tjdstart, double *dgeo, double *datm, double *dobs, int32 TypeEvent, int32 helflag, double *dret, char *serr_ret)
{
  double tjd, trise;
  char serr[AS_MAXCH];
  char ObjectName[30];
  int32 ipl, retval, helflag2, direct;
  dret[0] = tjdstart;
  if (TypeEvent == 1 || TypeEvent == 2) {
    if (serr_ret != NULL)
      strcpy(serr_ret, "error: the moon has no morning first or evening last");
    return ERR;
  }
  strcpy(ObjectName, "moon");
  ipl = SE_MOON;
  helflag2 = helflag;
  helflag2 &= ~SE_HELFLAG_HIGH_PRECISION;

  tjd = tjdstart - 30;

  if ((retval = find_conjunct_sun(tjd, ipl, helflag, TypeEvent, &tjd, serr)) == ERR)
    return ERR;

  retval = get_heliacal_day(tjd, dgeo, datm, dobs, ObjectName, helflag2, TypeEvent, &tjd, serr);
  if (retval != OK)
    goto moon_event_err;
  dret[0] = tjd;

  retval = time_optimum_visibility(tjd, dgeo, datm, dobs, ObjectName, helflag, &tjd, serr);
  if (retval == ERR) goto moon_event_err;
  dret[1] = tjd;

  direct = 1;
  if (TypeEvent == 4)
    direct = -1;
  retval = time_limit_invisible(tjd, dgeo, datm, dobs, ObjectName, helflag, direct, &tjd, serr);
  if (retval == ERR) goto moon_event_err;
  dret[2] = tjd;

  direct *= -1;
  retval = time_limit_invisible(dret[1], dgeo, datm, dobs, ObjectName, helflag, direct, &tjd, serr);
  dret[0] = tjd;
  if (retval == ERR) goto moon_event_err;
#if 1

  if (TypeEvent == 3) {
    if ((retval = my_rise_trans(tjd, SE_SUN, "", SE_CALC_SET, helflag, dgeo, datm, &trise, serr)) == ERR)
      return ERR;
    if (trise < dret[1]) {
      dret[0] = trise;

    }

  } else {
    if ((retval = my_rise_trans(dret[1], SE_SUN, "", SE_CALC_RISE, helflag, dgeo, datm, &trise, serr)) == ERR)
      return ERR;
    if (dret[0] > trise) {
      dret[0] = trise;

    }
  }
#endif

  if (TypeEvent == 4) {
    tjd = dret[0];
    dret[0] = dret[2];
    dret[2] = tjd;
  }
moon_event_err:
  if (serr_ret != NULL && *serr != '\0')
    strcpy(serr_ret, serr);
  return retval;
}

static int32 MoonEventJDut(double JDNDaysUTStart, double *dgeo, double *datm, double *dobs, int32 TypeEvent, int32 helflag, double *dret, char *serr)
{
  int32 avkind = helflag & SE_HELFLAG_AVKIND;
  if (avkind)
    return moon_event_arc_vis(JDNDaysUTStart, dgeo, datm, dobs, TypeEvent, helflag, dret, serr);
  else
    return moon_event_vis_lim(JDNDaysUTStart, dgeo, datm, dobs, TypeEvent, helflag, dret, serr);
}

static int32 heliacal_ut(double JDNDaysUTStart, double *dgeo, double *datm, double *dobs, char *ObjectName, int32 TypeEventIn, int32 helflag, double *dret, char *serr_ret)
{
  int32 avkind = helflag & SE_HELFLAG_AVKIND;
  if (avkind)
    return heliacal_ut_arc_vis(JDNDaysUTStart, dgeo, datm, dobs, ObjectName, TypeEventIn, helflag, dret, serr_ret);
  else
    return heliacal_ut_vis_lim(JDNDaysUTStart, dgeo, datm, dobs, ObjectName, TypeEventIn, helflag, dret, serr_ret);
}

int32 CALL_CONV swe_heliacal_ut(double JDNDaysUTStart, double *dgeo, double *datm, double *dobs, char *ObjectNameIn, int32 TypeEvent, int32 helflag, double *dret, char *serr_ret)
{
  int32 retval, Planet;
  char ObjectName[AS_MAXCH], serr[AS_MAXCH], s[AS_MAXCH];
  double tjd0 = JDNDaysUTStart, tjd, dsynperiod, tjdmax, tadd;
  int32 MaxCountSynodicPeriod = MAX_COUNT_SYNPER;
  char *sevent[7] = {"", "morning first", "evening last", "evening first", "morning last", "acronychal rising", "acronychal setting"};
  if (dgeo[2] < SEI_ECL_GEOALT_MIN || dgeo[2] > SEI_ECL_GEOALT_MAX) {
    if (serr_ret != NULL)
      sprintf(serr_ret, "location for heliacal events must be between %.0f and %.0f m above sea\n", SEI_ECL_GEOALT_MIN, SEI_ECL_GEOALT_MAX);
    return ERR;
  }
  swi_set_tid_acc(JDNDaysUTStart, helflag, 0, serr);
  if (helflag & SE_HELFLAG_LONG_SEARCH)
    MaxCountSynodicPeriod = MAX_COUNT_SYNPER_MAX;

  if (serr_ret != NULL)
    *serr_ret = '\0';

  strcpy_VBsafe(ObjectName, ObjectNameIn);
  tolower_string_star(ObjectName);
  default_heliacal_parameters(datm, dgeo, dobs, helflag);
  swe_set_topo(dgeo[0], dgeo[1], dgeo[2]);
  Planet = DeterObject(ObjectName);
  if (Planet == SE_SUN) {
    if (serr_ret != NULL) {
      strcpy(serr_ret, "the sun has no heliacal rising or setting\n");
    }
    return ERR;
  }

  if (Planet == SE_MOON) {
    if (TypeEvent == 1 || TypeEvent == 2) {
      if (serr_ret != NULL) {
        sprintf(serr_ret, "%s (event type %d) does not exist for the moon\n", sevent[TypeEvent], TypeEvent);
      }
      return ERR;
    }
    tjd = tjd0;
    retval = MoonEventJDut(tjd, dgeo, datm, dobs, TypeEvent, helflag, dret, serr);
    while (retval != -2 && *dret < tjd0) {
      tjd += 15;
      *serr = '\0';
      retval = MoonEventJDut(tjd, dgeo, datm, dobs, TypeEvent, helflag, dret, serr);
    }
    if (serr_ret != NULL && *serr != '\0')
      strcpy(serr_ret, serr);
    return retval;
  }

  if (!(helflag & SE_HELFLAG_AVKIND)) {
    if (Planet == -1 || Planet >= SE_MARS) {
      if (TypeEvent == 3 || TypeEvent == 4) {
	if (serr_ret != NULL) {
	  if (Planet == -1)
	    strcpy(s, ObjectName);
	  else
	    swe_get_planet_name(Planet, s);
	  sprintf(serr_ret, "%s (event type %d) does not exist for %s\n", sevent[TypeEvent], TypeEvent, s);
	}
	return ERR;
      }
    }
  }

  if (helflag & SE_HELFLAG_AVKIND) {
    if (Planet == -1 || Planet >= SE_MARS) {
      if (TypeEvent == SE_ACRONYCHAL_RISING)
	TypeEvent = 3;
      if (TypeEvent == SE_ACRONYCHAL_SETTING)
	TypeEvent = 4;
    }

  } else if (1) {
    if (TypeEvent == SE_ACRONYCHAL_RISING || TypeEvent == SE_ACRONYCHAL_SETTING) {
      if (serr_ret != NULL) {
	if (Planet == -1)
	  strcpy(s, ObjectName);
	else
	  swe_get_planet_name(Planet, s);
	sprintf(serr_ret, "%s (event type %d) is not provided for %s\n", sevent[TypeEvent], TypeEvent, s);
      }
      return ERR;
    }
  }
  dsynperiod = get_synodic_period(Planet);
  tjdmax = tjd0 + dsynperiod * MaxCountSynodicPeriod;
  tadd = dsynperiod * 0.6;
  if (Planet == SE_MERCURY)
    tadd = 30;

  retval = -2;
  for (tjd = tjd0; tjd < tjdmax && retval == -2; tjd += tadd) {
    *serr = '\0';
    retval = heliacal_ut(tjd, dgeo, datm, dobs, ObjectName, TypeEvent, helflag, dret, serr);

    while (retval != -2 && *dret < tjd0) {
      tjd += tadd;
      *serr = '\0';
      retval = heliacal_ut(tjd, dgeo, datm, dobs, ObjectName, TypeEvent, helflag, dret, serr);
    }
  }

  if ((helflag & SE_HELFLAG_SEARCH_1_PERIOD) && (retval == -2 || dret[0] > tjd0 + dsynperiod * 1.5)) {
    strcpy(serr, "no heliacal date found within this synodic period");
    retval = -2;
  } else if (retval == -2) {
    sprintf(serr, "no heliacal date found within %d synodic periods", MaxCountSynodicPeriod);
    retval = ERR;
  }
  if (serr_ret != NULL && *serr != '\0')
    strcpy(serr_ret, serr);
  return retval;
}
