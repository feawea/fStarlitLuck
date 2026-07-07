# include "swephexp.h"
# include "sweph.h"

static TLS AS_BOOL init_leapseconds_done = FALSE;

int CALL_CONV swe_date_conversion(int y,
		     int m,
		     int d,
		     double uttime,
		     char c,
		     double *tjd)
{
  int rday, rmon, ryear;
  double rut, jd;
  int gregflag = SE_JUL_CAL;
  if (c == 'g')
    gregflag = SE_GREG_CAL;
  rut = uttime;
  jd = swe_julday(y, m, d, rut, gregflag);
  swe_revjul(jd, gregflag, &ryear, &rmon, &rday, &rut);
  *tjd = jd;
  if (rmon == m && rday == d && ryear == y) {
    return OK;
  } else {
    return ERR;
  }
}

double CALL_CONV swe_julday(int year, int month, int day, double hour, int gregflag)
{
  double jd;
  double u,u0,u1,u2;
  u = year;
  if (month < 3) u -=1;
  u0 = u + 4712.0;
  u1 = month + 1.0;
  if (u1 < 4) u1 += 12.0;
  jd = floor(u0*365.25)
     + floor(30.6*u1+0.000001)
     + day + hour/24.0 - 63.5;
  if (gregflag == SE_GREG_CAL) {
    u2 = floor(fabs(u) / 100) - floor(fabs(u) / 400);
    if (u < 0.0) u2 = -u2;
    jd = jd - u2 + 2;
    if ((u < 0.0) && (u/100 == floor(u/100)) && (u/400 != floor(u/400)))
      jd -=1;
  }
  return jd;
}

void CALL_CONV swe_revjul (double jd, int gregflag,
	     int *jyear, int *jmon, int *jday, double *jut)
{
  double u0,u1,u2,u3,u4;
  u0 = jd + 32082.5;
  if (gregflag == SE_GREG_CAL) {
    u1 = u0 + floor (u0/36525.0) - floor (u0/146100.0) - 38.0;
    if (jd >= 1830691.5) u1 +=1;
    u0 = u0 + floor (u1/36525.0) - floor (u1/146100.0) - 38.0;
  }
  u2 = floor (u0 + 123.0);
  u3 = floor ( (u2 - 122.2) / 365.25);
  u4 = floor ( (u2 - floor (365.25 * u3) ) / 30.6001);
  *jmon = (int) (u4 - 1.0);
  if (*jmon > 12) *jmon -= 12;
  *jday = (int) (u2 - floor (365.25 * u3) - floor (30.6001 * u4));
  *jyear = (int) (u3 + floor ( (u4 - 2.0) / 12.0) - 4800);
  *jut = (jd - floor (jd + 0.5) + 0.5) * 24.0;
}

void CALL_CONV swe_utc_time_zone(
        int32 iyear, int32 imonth, int32 iday,
        int32 ihour, int32 imin, double dsec,
        double d_timezone,
        int32 *iyear_out, int32 *imonth_out, int32 *iday_out,
        int32 *ihour_out, int32 *imin_out, double *dsec_out
        )
{
  double tjd, d;
  AS_BOOL have_leapsec = FALSE;
  double dhour;
  if (dsec >= 60.0) {
    have_leapsec = TRUE;
    dsec -= 1.0;
  }
  dhour = ((double) ihour) + ((double) imin) / 60.0 + dsec / 3600.0;
  tjd = swe_julday(iyear, imonth, iday, 0, SE_GREG_CAL);
  dhour -= d_timezone;
  if (dhour < 0.0) {
    tjd -= 1.0;
    dhour += 24.0;
  }
  if (dhour >= 24.0) {
    tjd += 1.0;
    dhour -= 24.0;
  }
  swe_revjul(tjd + 0.001, SE_GREG_CAL, iyear_out, imonth_out, iday_out, &d);
  *ihour_out = (int) dhour;
  d = (dhour - (double) *ihour_out) * 60;
  *imin_out = (int) d;
  *dsec_out = (d - (double) *imin_out) * 60;
  if (have_leapsec)
    *dsec_out += 1.0;
}

#define NLEAP_SECONDS 27
#define NLEAP_SECONDS_SPACE 100
static TLS int leap_seconds[NLEAP_SECONDS_SPACE] = {
19720630,
19721231,
19731231,
19741231,
19751231,
19761231,
19771231,
19781231,
19791231,
19810630,
19820630,
19830630,
19850630,
19871231,
19891231,
19901231,
19920630,
19930630,
19940630,
19951231,
19970630,
19981231,
20051231,
20081231,
20120630,
20150630,
20161231,
0
};
#define J1972 2441317.5
#define NLEAP_INIT 10

static int init_leapsec(void)
{
  FILE *fp;
  int ndat, ndat_last;
  int tabsiz = 0;
  int i;
  char s[AS_MAXCH];
  char *sp;
  if (!init_leapseconds_done) {
    init_leapseconds_done = TRUE;
    tabsiz = NLEAP_SECONDS;
    ndat_last = leap_seconds[NLEAP_SECONDS - 1];

    if ((fp = swi_fopen(-1, "seleapsec.txt", swed.ephepath, NULL)) == NULL)
      return NLEAP_SECONDS;
    while(fgets(s, AS_MAXCH, fp) != NULL) {
      sp = s;
      while (*sp == ' ' || *sp == '\t') sp++;
        sp++;
      if (*sp == '#' || *sp == '\n')
        continue;
      ndat = atoi(s);
      if (ndat <= ndat_last)
        continue;

      if (tabsiz >= NLEAP_SECONDS_SPACE)
        return tabsiz;
      leap_seconds[tabsiz] = ndat;
      tabsiz++;
    }
    if (tabsiz > NLEAP_SECONDS) leap_seconds[tabsiz] = 0;
    fclose(fp);
    return tabsiz;
  }

  tabsiz = 0;
  for (i = 0; i < NLEAP_SECONDS_SPACE; i++) {
    if (leap_seconds[i] == 0)
      break;
    else
      tabsiz++;
  }
  return tabsiz;
}

int32 CALL_CONV swe_utc_to_jd(int32 iyear, int32 imonth, int32 iday, int32 ihour, int32 imin, double dsec, int32 gregflag, double *dret, char *serr)
{
  double tjd_ut1, tjd_et, tjd_et_1972, dhour, d;
  int iyear2, imonth2, iday2;
  int i, j, ndat, nleap, tabsiz_nleap;

  tjd_ut1 = swe_julday(iyear, imonth, iday, 0, gregflag);
  swe_revjul(tjd_ut1, gregflag, &iyear2, &imonth2, &iday2, &d);
  if (iyear != iyear2 || imonth != imonth2 || iday != iday2) {
    if (serr != NULL)
      sprintf(serr, "invalid date: year = %d, month = %d, day = %d", iyear, imonth, iday);
    return ERR;
  }
  if (ihour < 0 || ihour > 23
   || imin < 0 || imin > 59
   || dsec < 0 || dsec >= 61
   || (dsec >= 60 && (imin < 59 || ihour < 23 || tjd_ut1 < J1972))) {
    if (serr != NULL)
      sprintf(serr, "invalid time: %d:%d:%.2f", ihour, imin, dsec);
    return ERR;
  }
  dhour = (double) ihour + ((double) imin) / 60.0 + dsec / 3600.0;

  if (tjd_ut1 < J1972) {
    dret[1] = swe_julday(iyear, imonth, iday, dhour, gregflag);
    dret[0] = dret[1] + swe_deltat_ex(dret[1], -1, NULL);
    return OK;
  }

  if (gregflag == SE_JUL_CAL) {
    gregflag = SE_GREG_CAL;
    swe_revjul(tjd_ut1, gregflag, &iyear, &imonth, &iday, &d);
  }

  tabsiz_nleap = init_leapsec();
  nleap = NLEAP_INIT;
  ndat = iyear * 10000 + imonth * 100 + iday;
  for (i = 0; i < tabsiz_nleap; i++) {
    if (ndat <= leap_seconds[i])
      break;
    nleap++;
  }

  d = swe_deltat_ex(tjd_ut1, -1, NULL) * 86400.0;
  if (d - (double) nleap - 32.184 >= 1.0) {
    dret[1] = tjd_ut1 + dhour / 24.0;
    dret[0] = dret[1] + swe_deltat_ex(dret[1], -1, NULL);
    return OK;
  }

  if (dsec >= 60) {
    j = 0;
    for (i = 0; i < tabsiz_nleap; i++) {
      if (ndat == leap_seconds[i]) {
	j = 1;
	break;
      }
    }
    if (j != 1) {
      if (serr != NULL)
	sprintf(serr, "invalid time (no leap second!): %d:%d:%.2f", ihour, imin, dsec);
      return ERR;
    }
  }

  d = tjd_ut1 - J1972;

  d += (double) ihour / 24.0 + (double) imin / 1440.0 + dsec / 86400.0;

  tjd_et_1972 = J1972 + (32.184 + NLEAP_INIT) / 86400.0;
  tjd_et = tjd_et_1972 + d + ((double) (nleap - NLEAP_INIT)) / 86400.0;
  d = swe_deltat_ex(tjd_et, -1, NULL);
  tjd_ut1 = tjd_et - swe_deltat_ex(tjd_et - d, -1, NULL);
  tjd_ut1 = tjd_et - swe_deltat_ex(tjd_ut1, -1, NULL);
  dret[0] = tjd_et;
  dret[1] = tjd_ut1;
  return OK;
}

void CALL_CONV swe_jdet_to_utc(double tjd_et, int32 gregflag, int32 *iyear, int32 *imonth, int32 *iday, int32 *ihour, int32 *imin, double *dsec)
{
  int i;
  int second_60 = 0;
  int iyear2, imonth2, iday2, nleap, ndat, tabsiz_nleap;
  double d, tjd, tjd_et_1972, tjd_ut, dret[10];

  tjd_et_1972 = J1972 + (32.184 + NLEAP_INIT) / 86400.0;
  d = swe_deltat_ex(tjd_et, -1, NULL);
  tjd_ut = tjd_et - swe_deltat_ex(tjd_et - d, -1, NULL);
  tjd_ut = tjd_et - swe_deltat_ex(tjd_ut, -1, NULL);
  if (tjd_et < tjd_et_1972) {
    swe_revjul(tjd_ut, gregflag, iyear, imonth, iday, &d);
    *ihour = (int32) d;
    d -= (double) *ihour;
    d *= 60;
    *imin = (int32) d;
    *dsec = (d - (double) *imin) * 60.0;
    return;
  }

  tabsiz_nleap = init_leapsec();
  swe_revjul(tjd_ut-1, SE_GREG_CAL, &iyear2, &imonth2, &iday2, &d);
  ndat = iyear2 * 10000 + imonth2 * 100 + iday2;
  nleap = 0;
  for (i = 0; i < tabsiz_nleap; i++) {
    if (ndat <= leap_seconds[i])
      break;
    nleap++;
  }

  if (nleap < tabsiz_nleap) {
    i = leap_seconds[nleap];
    iyear2 = i / 10000;
    imonth2 = (i % 10000) / 100;;
    iday2 = i % 100;
    tjd = swe_julday(iyear2, imonth2, iday2, 0, SE_GREG_CAL);
    swe_revjul(tjd+1, SE_GREG_CAL, &iyear2, &imonth2, &iday2, &d);
    swe_utc_to_jd(iyear2,imonth2,iday2, 0, 0, 0, SE_GREG_CAL, dret, NULL);
    d = tjd_et - dret[0];
    if (d >= 0) {
      nleap++;
    } else if (d < 0 && d > -1.0/86400.0) {
      second_60 = 1;
    }
  }

  tjd = J1972 + (tjd_et - tjd_et_1972) - ((double) nleap + second_60) / 86400.0;
  swe_revjul(tjd, SE_GREG_CAL, iyear, imonth, iday, &d);
  *ihour = (int32) d;
  d -= (double) *ihour;
  d *= 60;
  *imin = (int32) d;
  *dsec = (d - (double) *imin) * 60.0 + second_60;

  d = swe_deltat_ex(tjd_et, -1, NULL);
  d = swe_deltat_ex(tjd_et - d, -1, NULL);
  if (d * 86400.0 - (double) (nleap + NLEAP_INIT) - 32.184 >= 1.0) {
    swe_revjul(tjd_et - d, SE_GREG_CAL, iyear, imonth, iday, &d);
    *ihour = (int32) d;
    d -= (double) *ihour;
    d *= 60;
    *imin = (int32) d;
    *dsec = (d - (double) *imin) * 60.0;
  }
  if (gregflag == SE_JUL_CAL) {
    tjd = swe_julday(*iyear, *imonth, *iday, 0, SE_GREG_CAL);
    swe_revjul(tjd, gregflag, iyear, imonth, iday, &d);
  }
}

void CALL_CONV swe_jdut1_to_utc(double tjd_ut, int32 gregflag, int32 *iyear, int32 *imonth, int32 *iday, int32 *ihour, int32 *imin, double *dsec)
{
  double tjd_et = tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL);
  swe_jdet_to_utc(tjd_et, gregflag, iyear, imonth, iday, ihour, imin, dsec);
}
