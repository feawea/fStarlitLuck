#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"
#include "swehouse.h"
#include <string.h>

#define MILLIARCSEC 	(1.0 / 3600000.0)
#define SOLAR_YEAR   365.24219893
#define ARMCS ((SOLAR_YEAR+1) / SOLAR_YEAR * 360)

static double Asc1(double, double, double, double);
static double AscDash(double, double, double, double);
static double Asc2(double, double, double, double);
static int CalcH(double th, double fi, double ekl, char hsy, struct houses *hsp);
static int sidereal_houses_ecl_t0(double tjde,
                           double armc,
                           double eps,
                           double *nutlo,
                           double lat,
			   int hsys,
                           double *cusp,
                           double *ascmc,
			   double *cusp_speed,
			   double *ascmc_speed,
			   char *serr);
static int sidereal_houses_trad(double tjde,
			   int32 iflag,
                           double armc,
                           double eps,
                           double nutl,
                           double lat,
			   int hsys,
                           double *cusp,
                           double *ascmc,
			   double *cusp_speed,
			   double *ascmc_speed,
			   char *serr);
static int sidereal_houses_ssypl(double tjde,
                           double armc,
                           double eps,
                           double *nutlo,
                           double lat,
			   int hsys,
                           double *cusp,
                           double *ascmc,
			   double *cusp_speed,
			   double *ascmc_speed,
			   char *serr);
static int sunshine_solution_makransky(double ramc, double lat, double ecl, struct houses *hsp);
static int sunshine_solution_treindl(double ramc, double lat, double ecl, struct houses *hsp);
#if 0
static void test_Asc1();
#endif

int CALL_CONV swe_houses(double tjd_ut,
				double geolat,
				double geolon,
				int hsys,
				double *cusp,
				double *ascmc)
{
  int i, retc = 0;
  double armc, eps, nutlo[2];
  double tjde = tjd_ut + swe_deltat_ex(tjd_ut, -1, NULL);
  eps = swi_epsiln(tjde, 0) * RADTODEG;
  swi_nutation(tjde, 0, nutlo);
  for (i = 0; i < 2; i++)
    nutlo[i] *= RADTODEG;
  armc = swe_degnorm(swe_sidtime0(tjd_ut, eps + nutlo[1], nutlo[0]) * 15 + geolon);
  if (toupper(hsys) ==  'I') {
    int flags = SEFLG_SPEED| SEFLG_EQUATORIAL;
    double xp[6];
    int result = swe_calc_ut(tjd_ut, SE_SUN, flags, xp, NULL);
    if (result < 0) {

      result = swe_houses_armc_ex2(armc, geolat, eps + nutlo[1], 'O', cusp, ascmc, NULL, NULL, NULL);
      return ERR;
    }
    ascmc[9] = xp[1];
  }
#ifdef TRACE
  swi_open_trace(NULL);
  if (swi_trace_count <= TRACE_COUNT_MAX) {
    if (swi_fp_trace_c != NULL) {
      fputs("\n[SWE_HOUSES]\n", swi_fp_trace_c);
      fprintf(swi_fp_trace_c, "#if 0\n");
      fprintf(swi_fp_trace_c, "  tjd = %.9f;", tjd_ut);
      fprintf(swi_fp_trace_c, " geolon = %.9f;", geolon);
      fprintf(swi_fp_trace_c, " geolat = %.9f;", geolat);
      fprintf(swi_fp_trace_c, " hsys = %d;\n", hsys);
      fprintf(swi_fp_trace_c, "  retc = swe_houses(tjd, geolat, geolon, hsys, cusp, ascmc);\n");
      fprintf(swi_fp_trace_c, "  [ swe_houses calls swe_houses_armc as follows: ]\n");
      fprintf(swi_fp_trace_c, "#endif\n");
      fflush(swi_fp_trace_c);
    }
  }
#endif
  retc = swe_houses_armc_ex2(armc, geolat, eps + nutlo[1], hsys, cusp, ascmc, NULL, NULL, NULL);
  return retc;
}

int CALL_CONV swe_houses_ex(double tjd_ut,
                                int32 iflag,
				double geolat,
				double geolon,
				int hsys,
				double *cusp,
				double *ascmc)
{
  return swe_houses_ex2(tjd_ut, iflag, geolat, geolon, hsys, cusp, ascmc, NULL, NULL, NULL);
}

int CALL_CONV swe_houses_ex2(double tjd_ut,
                                int32 iflag,
				double geolat,
				double geolon,
				int hsys,
				double *cusp,
				double *ascmc,
			        double *cusp_speed,
				double *ascmc_speed,
				char *serr)
{
  int i, retc = 0;
  double armc, eps_mean, nutlo[2];
  double tjde = tjd_ut + swe_deltat_ex(tjd_ut, iflag, NULL);
  struct sid_data *sip = &swed.sidd;
  double xp[6];
  int retc_makr = 0;
  int ito;
  if (toupper(hsys) == 'G')
    ito = 36;
  else
    ito = 12;
  if ((iflag & SEFLG_SIDEREAL) && !swed.ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);
  eps_mean = swi_epsiln(tjde, 0) * RADTODEG;
  swi_nutation(tjde, 0, nutlo);
  for (i = 0; i < 2; i++)
    nutlo[i] *= RADTODEG;
  if (iflag & SEFLG_NONUT) {
    for (i = 0; i < 2; i++)
      nutlo[i] = 0;
  }
#ifdef TRACE
  swi_open_trace(NULL);
  if (swi_trace_count <= TRACE_COUNT_MAX) {
    if (swi_fp_trace_c != NULL) {
      fputs("\n[SWE_HOUSES_EX]\n", swi_fp_trace_c);
      fprintf(swi_fp_trace_c, "#if 0\n");
      fprintf(swi_fp_trace_c, "  tjd = %.9f;", tjd_ut);
      fprintf(swi_fp_trace_c, " iflag = %d;\n", iflag);
      fprintf(swi_fp_trace_c, " geolon = %.9f;", geolon);
      fprintf(swi_fp_trace_c, " geolat = %.9f;", geolat);
      fprintf(swi_fp_trace_c, " hsys = %d;\n", hsys);
      fprintf(swi_fp_trace_c, "  retc = swe_houses_ex(tjd, iflag, geolat, geolon, hsys, cusp, ascmc);\n");
      fprintf(swi_fp_trace_c, "  [ swe_houses calls swe_houses_armc as follows: ]\n");
      fprintf(swi_fp_trace_c, "#endif\n");
      fflush(swi_fp_trace_c);
    }
  }
#endif

  armc = swe_degnorm(swe_sidtime0(tjd_ut, eps_mean + nutlo[1], nutlo[0]) * 15 + geolon);

  if (toupper(hsys) ==  'I') {
    int flags = SEFLG_SPEED| SEFLG_EQUATORIAL;
    retc_makr = swe_calc_ut(tjd_ut, SE_SUN, flags, xp, NULL);
    if (retc_makr < 0) {

      hsys = (int) 'O';
    }
    ascmc[9] = xp[1];
  }
  if (iflag & SEFLG_SIDEREAL) {
    if (sip->sid_mode & SE_SIDBIT_ECL_T0)
      retc = sidereal_houses_ecl_t0(tjde, armc, eps_mean + nutlo[1], nutlo, geolat, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);
    else if (sip->sid_mode & SE_SIDBIT_SSY_PLANE)
      retc = sidereal_houses_ssypl(tjde, armc, eps_mean + nutlo[1], nutlo, geolat, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);
    else
      retc = sidereal_houses_trad(tjde, iflag, armc, eps_mean + nutlo[1], nutlo[0], geolat, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);
  } else {
    retc = swe_houses_armc_ex2(armc, geolat, eps_mean + nutlo[1], hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);
    if (toupper(hsys) ==  'I')
      ascmc[9] = xp[1];
  }
  if (iflag & SEFLG_RADIANS) {
    for (i = 1; i <= ito; i++)
      cusp[i] *= DEGTORAD;
    for (i = 0; i < SE_NASCMC; i++)
      ascmc[i] *= DEGTORAD;
  }
  if (retc_makr < 0)
    return retc_makr;
  return retc;
}

static int sidereal_houses_ecl_t0(double tjde,
                           double armc,
                           double eps,
                           double *nutlo,
                           double lat,
			   int hsys,
                           double *cusp,
                           double *ascmc,
			   double *cusp_speed,
			   double *ascmc_speed,
			   char *serr)
{
  int i, j, retc = OK;
  double x[6], xvpx[6], x2[6], epst0, xnorm[6];
  double rxy, rxyz, c2, epsx, sgn, fac, dvpx, dvpxe;
  double armcx;
  struct sid_data *sip = &swed.sidd;
  int ito;
  if (toupper(hsys) == 'G')
    ito = 36;
  else
    ito = 12;

  epst0 = swi_epsiln(sip->t0, 0);

  x[0] = x[4] = 1;
  x[1] = x[2] = x[3] = x[5] = 0;

  swi_coortrf(x, x, -epst0);
  swi_coortrf(x+3, x+3, -epst0);

  swi_precess(x, sip->t0, 0, J_TO_J2000);
  swi_precess(x, tjde, 0, J2000_TO_J);
  swi_precess(x+3, sip->t0, 0, J_TO_J2000);
  swi_precess(x+3, tjde, 0, J2000_TO_J);

  swi_coortrf(x, x, (eps - nutlo[1]) * DEGTORAD);
  swi_coortrf(x+3, x+3, (eps - nutlo[1]) * DEGTORAD);
  swi_cartpol_sp(x, x);
  x[0] += nutlo[0] * DEGTORAD;
  swi_polcart_sp(x, x);
  swi_coortrf(x, x, -eps * DEGTORAD);
  swi_coortrf(x+3, x+3, -eps * DEGTORAD);

  swi_cross_prod(x, x+3, xnorm);
  rxy =  xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
  c2 = (rxy + xnorm[2] * xnorm[2]);
  rxyz = sqrt(c2);
  rxy = sqrt(rxy);
  epsx = asin(rxy / rxyz) * RADTODEG;

  if (fabs(x[5]) < 1e-15)
    x[5] = 1e-15;
  fac = x[2] / x[5];
  sgn = x[5] / fabs(x[5]);
  for (j = 0; j <= 2; j++)
    xvpx[j] = (x[j] - fac * x[j+3]) * sgn;

  swi_cartpol(xvpx, x2);
  dvpx = x2[0] * RADTODEG;

  armcx = swe_degnorm(armc - dvpx);

  retc = swe_houses_armc_ex2(armcx, lat, epsx, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);

  dvpxe = acos(swi_dot_prod_unit(x, xvpx)) * RADTODEG;
  if (tjde < sip->t0)
    dvpxe = -dvpxe;
  for (i = 1; i <= ito; i++)
    cusp[i] = swe_degnorm(cusp[i] - dvpxe - sip->ayan_t0);
  for (i = 0; i <= SE_NASCMC; i++) {
    if (i == 2)
      continue;
    ascmc[i] = swe_degnorm(ascmc[i] - dvpxe - sip->ayan_t0);
  }
  if (hsys == 'N') {
    for (i = 1; i <= ito; i++) {
      cusp[i] = (i - 1) * 30;
    }
  }
  return retc;
}

static int sidereal_houses_ssypl(double tjde,
                           double armc,
                           double eps,
                           double *nutlo,
                           double lat,
			   int hsys,
                           double *cusp,
                           double *ascmc,
			   double *cusp_speed,
			   double *ascmc_speed,
			   char *serr)
{
  int i, j, retc = OK;
  double x[6], x0[6], xvpx[6], x2[6], xnorm[6];
  double rxy, rxyz, c2, epsx, eps2000, sgn, fac, dvpx, dvpxe, x00;
  double armcx;
  struct sid_data *sip = &swed.sidd;
  int ito;
  if (toupper(hsys) == 'G')
    ito = 36;
  else
    ito = 12;
  eps2000 = swi_epsiln(J2000, 0);

  x[0] = x[4] = 1;
  x[1] = x[2] = x[3] = x[5] = 0;

  swi_coortrf(x, x, -SSY_PLANE_INCL);
  swi_coortrf(x+3, x+3, -SSY_PLANE_INCL);
  swi_cartpol_sp(x, x);
  x[0] += SSY_PLANE_NODE_E2000;
  swi_polcart_sp(x, x);

  swi_coortrf(x, x, -eps2000);
  swi_coortrf(x+3, x+3, -eps2000);

  swi_precess(x, tjde, 0, J2000_TO_J);
  swi_precess(x+3, tjde, 0, J2000_TO_J);

  swi_coortrf(x, x, (eps - nutlo[1]) * DEGTORAD);
  swi_coortrf(x+3, x+3, (eps - nutlo[1]) * DEGTORAD);
  swi_cartpol_sp(x, x);
  x[0] += nutlo[0] * DEGTORAD;
  swi_polcart_sp(x, x);
  swi_coortrf(x, x, -eps * DEGTORAD);
  swi_coortrf(x+3, x+3, -eps * DEGTORAD);

  swi_cross_prod(x, x+3, xnorm);
  rxy =  xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
  c2 = (rxy + xnorm[2] * xnorm[2]);
  rxyz = sqrt(c2);
  rxy = sqrt(rxy);
  epsx = asin(rxy / rxyz) * RADTODEG;

  if (fabs(x[5]) < 1e-15)
    x[5] = 1e-15;
  fac = x[2] / x[5];
  sgn = x[5] / fabs(x[5]);
  for (j = 0; j <= 2; j++)
    xvpx[j] = (x[j] - fac * x[j+3]) * sgn;

  swi_cartpol(xvpx, x2);
  dvpx = x2[0] * RADTODEG;

  armcx = swe_degnorm(armc - dvpx);

  retc = swe_houses_armc_ex2(armcx, lat, epsx, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);

  dvpxe = acos(swi_dot_prod_unit(x, xvpx)) * RADTODEG;

  dvpxe -= SSY_PLANE_NODE * RADTODEG;

  x0[0] = 1;
  x0[1] = x0[2] = 0;

  if (sip->t0 != J2000)
    swi_precess(x0, sip->t0, 0, J_TO_J2000);

  swi_coortrf(x0, x0, eps2000);

  swi_cartpol(x0, x0);
  x0[0] -= SSY_PLANE_NODE_E2000;
  swi_polcart(x0, x0);
  swi_coortrf(x0, x0, SSY_PLANE_INCL);
  swi_cartpol(x0, x0);
  x0[0] += SSY_PLANE_NODE;
  x00 = x0[0] * RADTODEG;
  for (i = 1; i <= ito; i++)
    cusp[i] = swe_degnorm(cusp[i] - dvpxe - sip->ayan_t0 - x00);
  for (i = 0; i <= SE_NASCMC; i++) {
    if (i == 2)
      continue;
    ascmc[i] = swe_degnorm(ascmc[i] - dvpxe - sip->ayan_t0 - x00);
  }
  if (hsys == 'N') {
    for (i = 1; i <= ito; i++) {
      cusp[i] = (i - 1) * 30;
    }
  }
  return retc;
}

static int sidereal_houses_trad(double tjde,
			   int32 iflag,
                           double armc,
                           double eps,
                           double nutl,
                           double lat,
			   int hsys,
                           double *cusp,
                           double *ascmc,
			   double *cusp_speed,
			   double *ascmc_speed,
			   char *serr)
{
  int i, retc = OK;
  double ay;
  int ito;
  int ihs = toupper(hsys);
  int ihs2 = ihs;

  retc = swe_get_ayanamsa_ex(tjde, iflag, &ay, NULL);

  if (ihs == 'G')
    ito = 36;
  else
    ito = 12;
  if (ihs == 'W')
    ihs2 = 'E';

  retc = swe_houses_armc_ex2(armc, lat, eps, ihs2, cusp, ascmc, cusp_speed, ascmc_speed, serr);

  for (i = 1; i <= ito; i++) {

    cusp[i] = swe_degnorm(cusp[i] - ay);
    if (ihs == 'W')
      cusp[i] -= fmod(cusp[i], 30);
  }
  if (ihs == 'N') {
    for (i = 1; i <= ito; i++) {
      cusp[i] = (i - 1) * 30;
    }
  }
  for (i = 0; i < SE_NASCMC; i++) {
    if (i == 2)
      continue;

    ascmc[i] = swe_degnorm(ascmc[i] - ay);
  }

  return retc;
}

int CALL_CONV swe_houses_armc(
				double armc,
				double geolat,
				double eps,
				int hsys,
				double *cusp,
				double *ascmc)
{
  return swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc, NULL, NULL, NULL);
}

int CALL_CONV swe_houses_armc_ex2(
				double armc,
				double geolat,
				double eps,
				int hsys,
				double *cusp,
				double *ascmc,
				double *cusp_speed,
				double *ascmc_speed,
				char *serr)
{
  struct houses h, hm1, hp1;
  int i, retc = 0, rm1, rp1;
  int ito;
  static double saved_sundec = 99;
  if (toupper(hsys) == 'G')
    ito = 36;
  else
    ito = 12;
  armc = swe_degnorm(armc);
  h.do_speed = FALSE;
  h.do_hspeed = FALSE;
  if (ascmc_speed != NULL || cusp_speed != NULL)
    h.do_speed = TRUE;
  if (cusp_speed != NULL)
    h.do_hspeed = TRUE;
  if (toupper(hsys) ==  'I') {
    if (ascmc[9] == 99) {
      h.sundec = 0;
      if (saved_sundec != 99) h.sundec = saved_sundec;
    } else {
      h.sundec = ascmc[9];
      saved_sundec = h.sundec;
    }
    if (h.sundec < -24 || h.sundec > 24) {
      sprintf(serr, "House system I (Sunshine) needs valid Sun declination in ascmc[9]");
      return ERR;
    }
  }
  retc = CalcH(armc, geolat, eps, (char)hsys, &h);
  cusp[0] = 0;
  if (h.do_hspeed) cusp_speed[0] = 0;

  if (retc < 0) {
    ito = 12;
    if (serr != NULL) strcpy(serr, h.serr);
  }
  for (i = 1; i <= ito; i++) {
    cusp[i] = h.cusp[i];
    if (h.do_hspeed) cusp_speed[i] = h.cusp_speed[i];
  }
  ascmc[0] = h.ac;
  ascmc[1] = h.mc;
  ascmc[2] = armc;
  ascmc[3] = h.vertex;
  ascmc[4] = h.equasc;
  ascmc[5] = h.coasc1;
  ascmc[6] = h.coasc2;
  ascmc[7] = h.polasc;
  for (i = SE_NASCMC; i < 10; i++)
    ascmc[i] = 0;
  if (toupper(hsys) ==  'I')
    ascmc[9] = h.sundec ;
  if (h.do_speed && ascmc_speed != NULL) {
    ascmc_speed[0] = h.ac_speed;
    ascmc_speed[1] = h.mc_speed;
    ascmc_speed[2] = h.armc_speed;
    ascmc_speed[3] = h.vertex_speed;
    ascmc_speed[4] = h.equasc_speed;
    ascmc_speed[5] = h.coasc1_speed;
    ascmc_speed[6] = h.coasc2_speed;
    ascmc_speed[7] = h.polasc_speed;
    for (i = SE_NASCMC; i < 10; i++)
      ascmc_speed[i] = 0;
  }
  if (h.do_interpol) {
    double dt = 1.0 / 86400;
    double darmc = dt * ARMCS;
    hm1.do_speed = FALSE;
    hm1.do_hspeed = FALSE;
    hp1.do_speed = FALSE;
    hp1.do_hspeed = FALSE;
    if (toupper(hsys) ==  'I') {
      hm1.sundec = h.sundec;
      hp1.sundec = h.sundec;
    }
    rm1 = CalcH(armc - darmc, geolat, eps, (char)hsys, &hm1);
    rp1 = CalcH(armc + darmc, geolat, eps, (char)hsys, &hp1);
    if (rp1 >= 0 && rm1 >=0) {
      if (fabs(swe_difdeg2n(hp1.ac, h.ac)) > 90) {
	hp1 = h;
	dt = dt / 2;
      } else if (fabs(swe_difdeg2n(hm1.ac, h.ac)) > 90) {
	hm1 = h;
	dt = dt / 2;
      }
      for (i = 1; i <= 12; i++) {
	double dx = swe_difdeg2n(hp1.cusp[i], hm1.cusp[i]);
	cusp_speed[i] = dx / 2 / dt ;
      }
    }
  }
#ifdef TRACE
  swi_open_trace(NULL);
  if (swi_trace_count <= TRACE_COUNT_MAX) {
    if (swi_fp_trace_c != NULL) {
      fputs("\n[SWE_HOUSES_ARMC_EX2]\n", swi_fp_trace_c);
      fprintf(swi_fp_trace_c, "  armc = %.9f;", armc);
      fprintf(swi_fp_trace_c, " geolat = %.9f;", geolat);
      fprintf(swi_fp_trace_c, " eps = %.9f;", eps);
      fprintf(swi_fp_trace_c, " hsys = %d;\n", hsys);
      fprintf(swi_fp_trace_c, "  retc = swe_houses_armc_ex2(armc, geolat, eps, hsys, cusp, ascmc, cusp_speed, ascmc_speed, serr);\n");
      fputs("  printf(\"swe_houses_armc_ex2: %f\\t%f\\t%f\\t%c\\t\\n\", ", swi_fp_trace_c);
      fputs("  armc, geolat, eps, hsys);\n", swi_fp_trace_c);
      fputs("  printf(\"retc = %d\\n\", retc);\n", swi_fp_trace_c);
      fputs("  printf(\"cusp:\\n\");\n", swi_fp_trace_c);
      fputs("  for (i = 1; i <= 12; i++)\n", swi_fp_trace_c);
      fputs("    printf(\"  %d\\t%f\\n\", i, cusp[i]);\n", swi_fp_trace_c);
      fputs("  printf(\"ascmc:\\n\");\n", swi_fp_trace_c);
      fputs("  for (i = 0; i < 10; i++)\n", swi_fp_trace_c);
      fputs("    printf(\"  %d\\t%f\\n\", i, ascmc[i]);\n", swi_fp_trace_c);
      fputs("  printf(\"cusp_speed:\\n\");\n", swi_fp_trace_c);
      fputs("  for (i = 1; i <= 12; i++)\n", swi_fp_trace_c);
      fputs("    printf(\"  %d\\t%f\\n\", i, cusp_speed[i]);\n", swi_fp_trace_c);
      fputs("  printf(\"ascmc_speed:\\n\");\n", swi_fp_trace_c);
      fputs("  for (i = 0; i < 10; i++)\n", swi_fp_trace_c);
      fputs("    printf(\"  %d\\t%f\\n\", i, ascmc_speed[i]);\n", swi_fp_trace_c);
      fflush(swi_fp_trace_c);
    }
    if (swi_fp_trace_out != NULL) {
      fprintf(swi_fp_trace_out, "swe_houses_armc_ex2: %f\t%f\t%f\t%c\t\n", armc, geolat, eps, hsys);
      fprintf(swi_fp_trace_out, "retc = %d\n", retc);
      fputs("cusp:\n", swi_fp_trace_out);
      for (i = 1; i <= 12; i++)
	fprintf(swi_fp_trace_out, "  %d\t%f\n", i, cusp[i]);
      fputs("ascmc:\n", swi_fp_trace_out);
      for (i = 0; i < 10; i++)
	fprintf(swi_fp_trace_out, "  %d\t%f\n", i, ascmc[i]);
      fflush(swi_fp_trace_out);
    }
  }
#endif
#if 0

for (i = 1; i <=12; i++) {
  double x[6];
  x[0] = cusp[i]; x[1] = 0; x[2] = 1;
  cusp[i] = (swe_house_pos(armc, geolat, eps, hsys, x, NULL) - 1) * 30;
}
#endif
  return retc;
}

static double apc_sector(int n, double ph, double e, double az)
{
   int k, is_below_hor = 0;
   double kv, a, dasc, dret;

   if (fabs(ph * RADTODEG) > 90 - VERY_SMALL) {
     kv = 0;
     dasc = 0;
   } else {
     kv   = atan(tan(ph) * tan(e) * cos(az)/(1 + tan(ph) * tan(e) * sin(az)));
     if (fabs(ph * RADTODEG) < VERY_SMALL) {
       dasc = (90 - VERY_SMALL) * DEGTORAD;
       if (ph < 0)
         dasc = -dasc;
     } else {
       dasc = atan(sin(kv) / tan(ph));
     }
   }

   if (n < 8) {
     is_below_hor = 1;
     k = n - 1;
   } else {
     k = n - 13;
   }

   if (is_below_hor) {
     a = kv + az + PI/2 + k * (PI/2 - kv) / 3;
   } else {
     a = kv + az + PI/2 + k * (PI/2 + kv) / 3;
   }
   a = swe_radnorm(a);
   dret = atan2(tan(dasc) * tan(ph) * sin(az) + sin(a),
      cos(e) * (tan(dasc) * tan(ph) * cos(az) + cos(a)) + sin(e) * tan(ph) * sin(az - a));
   dret = swe_degnorm(dret * RADTODEG);
   return dret;
}

const char *CALL_CONV swe_house_name(int hsys)
{
  int h = hsys;
  if (h != 'i') h = toupper(h);
  switch (h) {
  case 'A': return "equal";
  case 'B': return "Alcabitius";
  case 'C': return "Campanus";
  case 'D': return "equal (MC)";
  case 'E': return "equal";
  case 'F': return "Carter poli-equ.";
  case 'G': return "Gauquelin sectors";
  case 'H': return "horizon/azimut";
  case 'I': return "Sunshine";
  case 'i': return "Sunshine/alt.";
  case 'J': return "Savard-A";
  case 'K': return "Koch";
  case 'L': return "Pullen SD";
  case 'M': return "Morinus";
  case 'N': return "equal/1=Aries";
  case 'O': return "Porphyry";
  case 'Q': return "Pullen SR";
  case 'R': return "Regiomontanus";
  case 'S': return "Sripati";
  case 'T': return "Polich/Page";
  case 'U': return "Krusinski-Pisa-Goelzer";
  case 'V': return "equal/Vehlow";
  case 'W': return "equal/ whole sign";
  case 'X': return "axial rotation system/Meridian houses";
  case 'Y': return "APC houses";
  default: return "Placidus";
  }
}

#define SUNSHINE_KEEP_MC_SOUTH	0

double swi_armc_to_mc(double armc, double eps)
{
  double tant, mc;
  if (fabs(armc - 90) > VERY_SMALL
      && fabs(armc - 270) > VERY_SMALL) {
    tant = tand(armc);
    mc = atand(tant / cosd(eps));
    if (armc > 90 && armc <= 270)
      mc = swe_degnorm(mc + 180);
  } else {
    if (fabs(armc - 90) <= VERY_SMALL)
      mc = 90;
    else
      mc = 270;
  }
  return mc;
}

#define VERY_SMALL_PLAC_ITER (1.0 / 360000.0 )
static int CalcH(
	double th, double fi, double ekl, char hsy, struct houses *hsp)

{
  double tane, tanfi, cosfi, sinfi, tant, sina, cosa, th2;
  double a, c, f, fh1, fh2, xh1, xh2, xs1, xs2, rectasc, ad3, acmc, vemc;
  int 	i, ih, ih2, retc = OK;
  double sine, cose;
  double x[3], krHorizonLon;
  int niter_max = 100;
  double cuspsv;
  *hsp->serr = '\0';
  hsp->do_interpol = 0;
  cose  = cosd(ekl);
  sine  = sind(ekl);
  tane  = tand(ekl);

  if (fabs(fabs(fi) - 90) < VERY_SMALL) {
    if (fi < 0)
      fi = -90 + VERY_SMALL;
    else
      fi = 90 - VERY_SMALL;
  }
  tanfi = tand(fi);

  if (fabs(th - 90) > VERY_SMALL
      && fabs(th - 270) > VERY_SMALL) {
    tant = tand(th);
    hsp->mc = atand(tant / cose);
    if (th > 90 && th <= 270)
      hsp->mc = swe_degnorm(hsp->mc + 180);
  } else {
    if (fabs(th - 90) <= VERY_SMALL)
      hsp->mc = 90;
    else
      hsp->mc = 270;
  }
  hsp->mc = swe_degnorm(hsp->mc);
  if (hsp->do_speed) hsp->mc_speed = AscDash(th, 0, sine, cose);

  hsp->ac = Asc1(th + 90, fi, sine, cose);
  if (hsp->do_speed)
    hsp->ac_speed = AscDash(th + 90, fi, sine, cose);
  if (hsp->do_hspeed) {
    for (i = 0; i <= 12; i++)
      hsp->cusp_speed[i] = 0;
  }
  hsp->armc_speed = ARMCS;

  hsp->cusp[1] = hsp->ac;
  hsp->cusp[10] = hsp->mc;
  if (hsp->do_hspeed) {
    hsp->cusp_speed[1] = hsp->ac_speed;
    hsp->cusp_speed[10] = hsp->mc_speed;
  }

  if (hsy > 95 && hsy != 'i') {
    sprintf(hsp->serr, "use of lower case letters like %c for house systems is deprecated", hsy);
    hsy = (char) (hsy - 32);
  }
  switch (hsy) {
  case 'A':
  case 'E':
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
      hsp->cusp[1] = hsp->ac;
    }
    for (i = 2; i <=12; i++) {
      hsp->cusp[i] = swe_degnorm(hsp->cusp[1] + (i-1) * 30);
    }
    if (hsp->do_hspeed) {
      for (i = 1; i <=12; i++) {
	hsp->cusp_speed[i] = hsp->ac_speed;
      }
    }
    break;
  case 'D':
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
    }
    hsp->cusp[10] = hsp->mc;
    for (i = 11; i <= 12; i++)
      hsp->cusp[i] = swe_degnorm(hsp->cusp[10] + (i-10) * 30);
    for (i = 1; i <= 9; i++)
      hsp->cusp[i] = swe_degnorm(hsp->cusp[10] + (i + 2) * 30);
    if (hsp->do_hspeed) {
      for (i = 1; i <=12; i++) {
	hsp->cusp_speed[i] = hsp->mc_speed;
      }
    }
    break;
  case 'C':

    fh1 = asind(sind (fi) / 2);

    fh2 = asind(sqrt (3.0) / 2 * sind(fi));
    cosfi = cosd(fi);
    if (fabs(cosfi) == 0) {
      if (fi > 0)
	xh1 = xh2 = 90;
      else
	xh1 = xh2 = 270;
    } else {

      xh1 = atand(sqrt (3.0) / cosfi);

      xh2 = atand(1 / sqrt (3.0) / cosfi);
    }
    hsp->cusp[11] = Asc1(th + 90 - xh1, fh1, sine, cose);
    hsp->cusp[12] = Asc1(th + 90 - xh2, fh2, sine, cose);
    hsp->cusp[2] = Asc1(th + 90 + xh2, fh2, sine, cose);
    hsp->cusp[3] = Asc1(th + 90 + xh1, fh1, sine, cose);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[11] = AscDash(th + 90 - xh1, fh1, sine, cose);
      hsp->cusp_speed[12] = AscDash(th + 90 - xh2, fh2, sine, cose);
      hsp->cusp_speed[2] = AscDash(th + 90 + xh2, fh2, sine, cose);
      hsp->cusp_speed[3] = AscDash(th + 90 + xh1, fh1, sine, cose);
    }

    if (fabs(fi) >= 90 - ekl) {
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
        hsp->ac = swe_degnorm(hsp->ac + 180);
        hsp->mc = swe_degnorm(hsp->mc + 180);
	for (i = 1; i <= 12; i++) {
	  if (i >= 4 && i < 10) continue;
	  hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
        }
      }
    }
    break;
  case 'H':
    if (fi > 0)
      fi = 90 - fi;
    else
      fi = -90 - fi;

    if (fabs(fabs(fi) - 90) < VERY_SMALL) {
      if (fi < 0)
	fi = -90 + VERY_SMALL;
      else
	fi = 90 - VERY_SMALL;
    }
    th = swe_degnorm(th + 180);
    fh1 = asind(sind (fi) / 2);
    fh2 = asind(sqrt (3.0) / 2 * sind(fi));
    cosfi = cosd(fi);
    if (fabs(cosfi) == 0) {
      if (fi > 0)
	xh1 = xh2 = 90;
      else
	xh1 = xh2 = 270;
    } else {

      xh1 = atand(sqrt (3.0) / cosfi);

      xh2 = atand(1 / sqrt (3.0) / cosfi);
    }
    hsp->cusp[11] = Asc1(th + 90 - xh1, fh1, sine, cose);
    hsp->cusp[12] = Asc1(th + 90 - xh2, fh2, sine, cose);
    hsp->cusp[1] = Asc1(th + 90, fi, sine, cose);
    hsp->cusp[2] = Asc1(th + 90 + xh2, fh2, sine, cose);
    hsp->cusp[3] = Asc1(th + 90 + xh1, fh1, sine, cose);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[11] = AscDash(th + 90 - xh1, fh1, sine, cose);
      hsp->cusp_speed[12] = AscDash(th + 90 - xh2, fh2, sine, cose);
      hsp->cusp_speed[1] = AscDash(th + 90, fi, sine, cose);
      hsp->cusp_speed[2] = AscDash(th + 90 + xh2, fh2, sine, cose);
      hsp->cusp_speed[3] = AscDash(th + 90 + xh1, fh1, sine, cose);
    }

    if (fabs(fi) >= 90 - ekl) {
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
        hsp->ac = swe_degnorm(hsp->ac + 180);
        hsp->mc = swe_degnorm(hsp->mc + 180);
	for (i = 1; i <= 12; i++) {
	  if (i >= 4 && i < 10) continue;
	  hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
        }
      }
    }
    for (i = 1; i <= 3; i++)
      hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
    for (i = 11; i <= 12; i++)
      hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);

    if (fi > 0)
      fi = 90 - fi;
    else
      fi = -90 - fi;
    th = swe_degnorm(th + 180);
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {
      hsp->ac = swe_degnorm(hsp->ac + 180);
    }
    break;
  case 'I':
  case 'i':
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
      hsp->cusp[1] = hsp->ac;
      if (! SUNSHINE_KEEP_MC_SOUTH && hsy == 'I') {
	hsp->mc = swe_degnorm(hsp->mc + 180);
	hsp->cusp[10] = hsp->mc;
      }
    }
    hsp->cusp[4] = swe_degnorm(hsp->cusp[10] + 180);
    hsp->cusp[7] = swe_degnorm(hsp->cusp[1] + 180);
    if (hsy == 'I') {
      retc = sunshine_solution_treindl(th, fi, ekl, hsp);
    } else {
      retc = sunshine_solution_makransky(th, fi, ekl, hsp);
    }
    if (retc == ERR) {
      strcpy(hsp->serr, "within polar circle, switched to Porphyry");
      hsy = 'O';
      goto porphyry;
    }
    hsp->do_interpol = hsp->do_hspeed;
    break;
  case 'J':

    sinfi = sind(fi);
    cosfi = cosd(fi);
    if (fabs(fi) < VERY_SMALL) {
      xs2 = 1 / 3.0;
      xs1 = 2 / 3.0;
    } else {
      xs2 = sind(fi / 3) / sinfi;
      xs1 = sind(2 * fi / 3) / sinfi;
    }
    xs2 = asind(xs2);
    xs1 = asind(xs1);

    if (cosfi == 0) {
      if (fi > 0)
	xh1 = xh2 = 90;
      else
	xh1 = xh2 = 270;
    } else {
      xh1 = atand(tand(xs1) / cosfi);
      xh2 = atand(tand(xs2) / cosfi);
    }

    fh1 = asind(sind (fi) * sind(90 - xs1));
    fh2 = asind(sind (fi) * sind(90 - xs2));
    hsp->cusp[12] = Asc1(th + 90 - xh2, fh2, sine, cose);
    hsp->cusp[11] = Asc1(th + 90 - xh1, fh1, sine, cose);
    hsp->cusp[2] = Asc1(th + 90 + xh2, fh2, sine, cose);
    hsp->cusp[3] = Asc1(th + 90 + xh1, fh1, sine, cose);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[11] = AscDash(th + 90 - xh1, fh1, sine, cose);
      hsp->cusp_speed[12] = AscDash(th + 90 - xh2, fh2, sine, cose);
      hsp->cusp_speed[3] = AscDash(th + 90 + xh1, fh1, sine, cose);
      hsp->cusp_speed[2] = AscDash(th + 90 + xh2, fh2, sine, cose);
    }

    if (fabs(fi) >= 90 - ekl) {
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
        hsp->ac = swe_degnorm(hsp->ac + 180);
        hsp->mc = swe_degnorm(hsp->mc + 180);
	for (i = 1; i <= 12; i++) {
	  if (i >= 4 && i < 10) continue;
	  hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
        }
      }
    }
    break;
  case 'K':
    if (fabs(fi) >= 90 - ekl) {
      retc = ERR;
      strcpy(hsp->serr, "within polar circle, switched to Porphyry");
      goto porphyry;
    }
    sina = sind(hsp->mc) * sine / cosd(fi);
    if (sina > 1) sina = 1;
    if (sina < -1) sina = -1;
    cosa = sqrt(1 - sina * sina);
    c = atand(tanfi / cosa);
    ad3 = asind(sind(c) * sina) / 3.0;
    hsp->cusp[11] = Asc1(th + 30 - 2 * ad3, fi, sine, cose);
    hsp->cusp[12] = Asc1(th + 60 - ad3, fi, sine, cose);
    hsp->cusp[2] = Asc1(th + 120 + ad3, fi, sine, cose);
    hsp->cusp[3] = Asc1(th + 150 + 2 * ad3, fi, sine, cose);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[11] = AscDash(th + 30 - 2 * ad3, fi, sine, cose);
      hsp->cusp_speed[12] = AscDash(th + 60 - ad3, fi, sine, cose);
      hsp->cusp_speed[2] = AscDash(th + 120 + ad3, fi, sine, cose);
      hsp->cusp_speed[3] = AscDash(th + 150 + 2 * ad3, fi, sine, cose);
    }
    break;
  case 'L':
    {
      double d, q1;
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {

	hsp->ac = swe_degnorm(hsp->ac + 180);
	hsp->cusp[1] = hsp->ac;
	acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      }
      q1 = 180 - acmc;
      d = (acmc - 90) / 4.0;
      if (acmc <= 30) {
	hsp->cusp[11] = hsp->cusp[12] = swe_degnorm(hsp->mc + acmc / 2);
      } else {
	hsp->cusp[11] = swe_degnorm(hsp->mc + 30 + d);
	hsp->cusp[12] = swe_degnorm(hsp->mc + 60 + 3 * d);
      }
      d = (q1 - 90) / 4.0;
      if (q1 <= 30) {
	hsp->cusp[2] = hsp->cusp[3] = swe_degnorm(hsp->ac + q1 / 2);
      } else {
	hsp->cusp[2] = swe_degnorm(hsp->ac + 30 + d);
	hsp->cusp[3] = swe_degnorm(hsp->ac + 60 + 3 * d);
      }
    }
    hsp->do_interpol = hsp->do_hspeed;
    break;
  case 'N':
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
    }
    for (i = 1; i <= 12; i++)
      hsp->cusp[i] = (i - 1) * 30.0;
    break;
  case 'O':
porphyry:
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
      hsp->cusp[1] = hsp->ac;
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    }
    hsp->cusp[1] = hsp->ac;
    hsp->cusp[10] = hsp->mc;
    hsp->cusp[2] = swe_degnorm(hsp->ac + (180 - acmc) / 3);
    hsp->cusp[3] = swe_degnorm(hsp->ac + (180 - acmc) / 3 * 2);
    hsp->cusp[11] = swe_degnorm(hsp->mc + acmc / 3);
    hsp->cusp[12] = swe_degnorm(hsp->mc + acmc / 3 * 2);
    if (hsp->do_hspeed) {
      double q1_speed = hsp->ac_speed - hsp->mc_speed;

      hsp->cusp_speed[1] = hsp->ac_speed;
      hsp->cusp_speed[10] = hsp->mc_speed;
      hsp->cusp_speed[2] = hsp->ac_speed  - q1_speed / 3;
      hsp->cusp_speed[3] = hsp->ac_speed  - q1_speed / 3 * 2;
      hsp->cusp_speed[11] = hsp->ac_speed  + q1_speed / 3;
      hsp->cusp_speed[12] = hsp->ac_speed  + q1_speed / 3 * 2;
    }
    break;
  case 'Q':
    {
      double q, c, csq, ccr, cqx, two23, third, r, r1, r2, x, xr, xr3, xr4;
      third = 1.0 / 3.0;
      two23 = pow(2.0 * 2.0, third);
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {

       hsp->ac = swe_degnorm(hsp->ac + 180);
       hsp->cusp[1] = hsp->ac;
       acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      }
      q = acmc;
      if (q > 90) q = 180 - q;
      if (q < 1e-30) {

	x = xr = xr3 = 0;
	xr4 = 180;
      } else {
	c = (180 - q) / q;
	csq = c * c;
	ccr = pow(csq - c, third);
	cqx = sqrt(two23 * ccr + 1.0);
	r1 = 0.5 * cqx;
	r2 = 0.5 * sqrt(-2*(1-2*c) / cqx - two23 * ccr + 2);
	r = r1 + r2 - 0.5;
	x = q / (2 * r + 1);
	xr = r * x;
	xr3 = xr * r * r;
	xr4 = xr3 * r;
      }
      if (acmc > 90) {
	hsp->cusp[11] = swe_degnorm(hsp->mc + xr3);
	hsp->cusp[12] = swe_degnorm(hsp->cusp[11] + xr4);
	hsp->cusp[2] = swe_degnorm(hsp->ac + xr);
	hsp->cusp[3] = swe_degnorm(hsp->cusp[2] + x);
      } else {
	hsp->cusp[11] = swe_degnorm(hsp->mc + xr);
	hsp->cusp[12] = swe_degnorm(hsp->cusp[11] + x);
	hsp->cusp[2] = swe_degnorm(hsp->ac + xr3);
	hsp->cusp[3] = swe_degnorm(hsp->cusp[2] + xr4);
      }
    }
    hsp->do_interpol = hsp->do_hspeed;
    break;
  case 'R':
    fh1 = atand (tanfi * 0.5);
    fh2 = atand (tanfi * cosd(30));
    hsp->cusp[11] = Asc1(30 + th, fh1, sine, cose);
    hsp->cusp[12] = Asc1(60 + th, fh2, sine, cose);
    hsp->cusp[2] = Asc1(120 + th, fh2, sine, cose);
    hsp->cusp[3] = Asc1(150 + th, fh1, sine, cose);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[11] = AscDash(30 + th, fh1, sine, cose);
      hsp->cusp_speed[12] = AscDash(60 + th, fh2, sine, cose);
      hsp->cusp_speed[2] = AscDash(120 + th, fh2, sine, cose);
      hsp->cusp_speed[3] = AscDash(150 + th, fh1, sine, cose);
    }

    if (fabs(fi) >= 90 - ekl) {
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
        hsp->ac = swe_degnorm(hsp->ac + 180);
        hsp->mc = swe_degnorm(hsp->mc + 180);
	for (i = 1; i <= 12; i++) {
	  if (i >= 4 && i < 10) continue;
	  hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
        }
      }
    }
    break;
  case 'S':

    {
      double s1, s4, q1;
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {

	hsp->ac = swe_degnorm(hsp->ac + 180);
	acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      }
      q1 = 180 - acmc;
      s1 = q1 / 3.0;
      s4 = acmc / 3.0;
      hsp->cusp[1] = swe_degnorm(hsp->ac - s4 * 0.5);
      hsp->cusp[2] = swe_degnorm(hsp->ac + s1 * 0.5);
      hsp->cusp[3] = swe_degnorm(hsp->ac + s1 * 1.5);
      hsp->cusp[10] = swe_degnorm(hsp->mc - s1 * 0.5);
      hsp->cusp[11] = swe_degnorm(hsp->mc + s4 * 0.5);
      hsp->cusp[12] = swe_degnorm(hsp->mc + s4 * 1.5);
    }
    hsp->do_interpol = hsp->do_hspeed;
    break;
  case 'T':
    fh1 = atand (tanfi / 3.0);
    fh2 = atand (tanfi * 2.0 / 3.0);
    hsp->cusp[11] =  Asc1(30 + th, fh1, sine, cose);
    hsp->cusp[12] =  Asc1(60 + th, fh2, sine, cose);
    hsp->cusp[2] =  Asc1(120 + th, fh2, sine, cose);
    hsp->cusp[3] =  Asc1(150 + th, fh1, sine, cose);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[11] =  AscDash(30 + th, fh1, sine, cose);
      hsp->cusp_speed[12] =  AscDash(60 + th, fh2, sine, cose);
      hsp->cusp_speed[2] =  AscDash(120 + th, fh2, sine, cose);
      hsp->cusp_speed[3] =  AscDash(150 + th, fh1, sine, cose);
    }

    if (fabs(fi) >= 90 - ekl) {
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
        hsp->ac = swe_degnorm(hsp->ac + 180);
	hsp->mc = swe_degnorm(hsp->mc + 180);
	for (i = 1; i <= 12; i++)
	  hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
      }
    }
    break;
  case 'V':
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
    }
    hsp->cusp[1] = swe_degnorm(hsp->ac - 15);
    for (i = 2; i <=12; i++)
      hsp->cusp[i] = swe_degnorm(hsp->cusp[1] + (i-1) * 30);
    if (hsp->do_hspeed) {
      for (i = 1; i <=12; i++) {
	hsp->cusp_speed[i] = hsp->ac_speed;
      }
    }
    break;
  case 'W':
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
      hsp->cusp[1] = hsp->ac;
    }
    hsp->cusp[1] = hsp->ac - fmod(hsp->ac, 30);
    for (i = 2; i <=12; i++)
      hsp->cusp[i] = swe_degnorm(hsp->cusp[1] + (i-1) * 30);
    break;
  case 'X': {

    int j;
    double a = th;
    for (i = 1; i <= 12; i++) {
      j = i + 10;
      if (j > 12) j -= 12;
      a = swe_degnorm(a + 30);
	  if (fabs(a - 90) > VERY_SMALL
        && fabs(a - 270) > VERY_SMALL) {
        tant = tand(a);
        hsp->cusp[j] = atand(tant / cose);
        if (a > 90 && a <= 270)
          hsp->cusp[j] = swe_degnorm(hsp->cusp[j] + 180);
      } else {
        if (fabs(a - 90) <= VERY_SMALL)
          hsp->cusp[j] = 90;
        else
		  hsp->cusp[j] = 270;
      }
	  hsp->cusp[j] = swe_degnorm(hsp->cusp[j]);
    }
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {
      hsp->ac = swe_degnorm(hsp->ac + 180);
    }
    hsp->do_interpol = hsp->do_hspeed;
    break; }
  case 'M': {

    int j;
    double a = th;
    double x[3];
    for (i = 1; i <= 12; i++) {
      j = i + 10;
      if (j > 12) j -= 12;
      a = swe_degnorm(a + 30);
      x[0] = a;
      x[1] = 0;
      swe_cotrans(x, x, ekl);
      hsp->cusp[j] = x[0];
    }
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {
      hsp->ac = swe_degnorm(hsp->ac + 180);
    }
    hsp->do_interpol = hsp->do_hspeed;
    break; }
  case 'F': {

    double a, ra;
    double x[3];
    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {

      hsp->ac = swe_degnorm(hsp->ac + 180);
      hsp->cusp[1] = hsp->ac;
    }
    x[0] = hsp->ac;
    x[1] = 0;
    swe_cotrans(x, x, -ekl);
    a = x[0];
    for (i = 2; i <= 12; i++) {
      if (i <= 3 || i >= 10) {
        ra = swe_degnorm(a + (i - 1) * 30);
	if (fabs(ra - 90) > VERY_SMALL
	  && fabs(ra - 270) > VERY_SMALL) {
	  tant = tand(ra);
	  hsp->cusp[i] = atand(tant / cose);
	  if (ra > 90 && ra <= 270)
	    hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
	} else {
	  if (fabs(ra - 90) <= VERY_SMALL)
	    hsp->cusp[i] = 90;
	  else
	    hsp->cusp[i] = 270;
	}
	hsp->cusp[i] = swe_degnorm(hsp->cusp[i]);
      }
    }
    hsp->do_interpol = hsp->do_hspeed;
    break; }
  case 'B': {

      double dek, r, sna, sda, sn3, sd3;
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
	hsp->ac = swe_degnorm(hsp->ac + 180);
	hsp->cusp[1] = hsp->ac;
	acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      }
      dek = asind(sind(hsp->ac) * sine);

      r = -tanfi * tand(dek);

      if (r > 1) r = 1;
      if (r < -1) r = -1;
      sda = acosd(r);
      sna = 180 - sda;
      sd3 = sda / 3;
      sn3 = sna / 3;
      rectasc = swe_degnorm(th + sd3);

      hsp->cusp[11] = Asc1(rectasc, 0, sine, cose);
      rectasc = swe_degnorm(th + 2 * sd3);
      hsp->cusp[12] = Asc1(rectasc, 0, sine, cose);
      rectasc = swe_degnorm(th + 180 - 2 * sn3);
      hsp->cusp[2] = Asc1(rectasc, 0, sine, cose);
      rectasc = swe_degnorm(th + 180 -  sn3);
      hsp->cusp[3] = Asc1(rectasc, 0, sine, cose);
    }
    hsp->do_interpol = hsp->do_hspeed;
    break;
  case 'G':
    for (i = 1; i <= 36; i++) {
      hsp->cusp[i] = 0;
      hsp->cusp_speed[i] = 0;
    }
    if (fabs(fi) >= 90 - ekl) {
      retc = ERR;
      strcpy(hsp->serr, "within polar circle, switched to Porphyry");
      hsy = (int) 'O';
      goto porphyry;
    }

    a = asind(tand(fi) * tane);
    for (ih = 2; ih <= 9; ih++) {
      ih2 = 10 - ih;
      fh1 = atand(sind(a * ih2 / 9) / tane);
      rectasc = swe_degnorm((90 / 9) * ih2 + th);
      tant = tand(asind(sine * sind(Asc1(rectasc, fh1, sine, cose))));
      if (fabs(tant) < VERY_SMALL) {
	hsp->cusp[ih] = rectasc;
	if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
      } else {

	f = atand(sind(asind(tanfi * tant) * ih2 / 9)  /tant);
        hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	cuspsv = 0;
        for (i = 1; i <= niter_max; i++) {
	  tant = tand(asind(sine * sind(hsp->cusp[ih])));
	  if (fabs(tant) < VERY_SMALL) {
	    hsp->cusp[ih] = rectasc;
	    if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
	    break;
	  }

	  f = atand(sind(asind(tanfi * tant) * ih2 / 9) / tant);
  	  hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	  if (i > 1 && fabs(swe_difdeg2n(hsp->cusp[ih], cuspsv)) < VERY_SMALL_PLAC_ITER)
	    break;
	  cuspsv = hsp->cusp[ih];
        }
#ifdef DEBUG_PLAC_ITER
  fprintf(stderr, "h=%d, niter=%d\n", ih, i);
#endif
	if (i >= niter_max) {
	  retc = ERR;
	  hsy = (int) 'O';
	  strcpy(hsp->serr, "very close to polar circle, switched to Porphyry");
	  goto porphyry;
	}
	if (hsp->do_hspeed) hsp->cusp_speed[ih] = AscDash(rectasc, f, sine, cose);
      }
      hsp->cusp[ih+18] = swe_degnorm(hsp->cusp[ih] + 180);
      if (hsp->do_hspeed) hsp->cusp_speed[ih + 18] = hsp->cusp_speed[ih];
    }

    for (ih = 29; ih <= 36; ih++) {
      ih2 = ih - 28;
      fh1 = atand(sind(a * ih2 / 9) / tane);
      rectasc = swe_degnorm(180 - ih2 * 90 / 9 + th);
      tant = tand(asind(sine * sind(Asc1(rectasc, fh1, sine, cose))));
      if (fabs(tant) < VERY_SMALL) {
        hsp->cusp[ih] = rectasc;
	if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
      } else {
        f = atand(sind(asind(tanfi * tant) * ih2 / 9) / tant);

        hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	cuspsv = 0;
        for (i = 1; i <= niter_max; i++) {
	  tant = tand(asind(sine * sind(hsp->cusp[ih])));
	  if (fabs(tant) < VERY_SMALL) {
	    hsp->cusp[ih] = rectasc;
	    if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
	    break;
	  }
	  f = atand(sind(asind(tanfi * tant) * ih2 / 9) / tant);

  	  hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	  if (i > 1 && fabs(swe_difdeg2n(hsp->cusp[ih], cuspsv)) < VERY_SMALL_PLAC_ITER)
	    break;
	  cuspsv = hsp->cusp[ih];
	}
#ifdef DEBUG_PLAC_ITER
  fprintf(stderr, "h=%d, niter=%d\n", ih, i);
#endif
	if (i >= niter_max) {
	  retc = ERR;
	  hsy = (int) 'O';
	  strcpy(hsp->serr, "very close to polar circle, switched to Porphyry");
	  goto porphyry;
	}
	if (hsp->do_hspeed) hsp->cusp_speed[ih] = AscDash(rectasc, f, sine, cose);
      }
      hsp->cusp[ih-18] = swe_degnorm(hsp->cusp[ih] + 180);
      if (hsp->do_hspeed) hsp->cusp_speed[ih - 18] = hsp->cusp_speed[ih];
    }
    hsp->cusp[1] = hsp->ac;
    hsp->cusp[10] = hsp->mc;
    hsp->cusp[19] = swe_degnorm(hsp->ac + 180);
    hsp->cusp[28] = swe_degnorm(hsp->mc + 180);
    if (hsp->do_hspeed) {
      hsp->cusp_speed[1] = hsp->ac_speed;
      hsp->cusp_speed[10] = hsp->mc_speed;
      hsp->cusp_speed[19] = hsp->ac_speed;
      hsp->cusp_speed[28] = hsp->mc_speed;
    }
    break;
  case 'U':

    acmc = swe_difdeg2n(hsp->ac, hsp->mc);
    if (acmc < 0) {
      hsp->ac = swe_degnorm(hsp->ac + 180);
    }

    x[0] = hsp->ac;
    x[1] = 0.0;
    x[2] = 1.0;
    swe_cotrans(x, x, -ekl);
    x[0] = x[0] - (th-90);
    swe_cotrans(x, x, -(90-fi));
    krHorizonLon = x[0];
    x[0] = x[0] - x[0];
    swe_cotrans(x, x, -90);

    for(i = 0; i < 6; i++) {

      x[0] = 30.0*i;
      x[1] = 0.0;
      swe_cotrans(x, x, 90);
      x[0] = x[0] + krHorizonLon;
      swe_cotrans(x, x, 90-fi);
      x[0] = swe_degnorm(x[0] + (th-90));

      hsp->cusp[i+1] = atand(tand(x[0])/cosd(ekl));
      if (x[0] > 90 && x[0] <= 270)
	hsp->cusp[i+1] = swe_degnorm(hsp->cusp[i+1] + 180);
      hsp->cusp[i+1] = swe_degnorm(hsp->cusp[i+1]);
      hsp->cusp[i+7] = swe_degnorm(hsp->cusp[i+1]+180);
    }
    break;
  case 'Y':
    for (i = 1; i <= 12; i++) {
      hsp->cusp[i] = apc_sector(i, fi * DEGTORAD, ekl * DEGTORAD, th * DEGTORAD);
    }

    hsp->cusp[10] = hsp->mc;
    hsp->cusp[4] = swe_degnorm(hsp->mc + 180);

    if (fabs(fi) >= 90 - ekl) {
      acmc = swe_difdeg2n(hsp->ac, hsp->mc);
      if (acmc < 0) {
        hsp->ac = swe_degnorm(hsp->ac + 180);
        hsp->mc = swe_degnorm(hsp->mc + 180);
	for (i = 1; i <= 12; i++)
	  hsp->cusp[i] = swe_degnorm(hsp->cusp[i] + 180);
      }
    }
    hsp->do_interpol = hsp->do_hspeed;
    break;
  default:
    if (fabs(fi) >= 90 - ekl) {
      retc = ERR;
      strcpy(hsp->serr, "within polar circle, switched to Porphyry");
      goto porphyry;
    }
    a = asind(tand(fi) * tane);
    fh1 = atand(sind(a / 3) / tane);
    fh2 = atand(sind(a * 2 / 3) / tane);

    rectasc = swe_degnorm(30 + th);
    tant = tand(asind(sine * sind(Asc1(rectasc, fh1, sine, cose))));
    ih = 11;
    if (fabs(tant) < VERY_SMALL) {
      hsp->cusp[ih] = rectasc;
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
    } else {

      f = atand(sind(asind(tanfi * tant) / 3)  /tant);
      hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
      cuspsv = 0;
      for (i = 1; i <= niter_max; i++) {
	tant = tand(asind(sine * sind(hsp->cusp[ih])));
	if (fabs(tant) < VERY_SMALL) {
	  hsp->cusp[ih] = rectasc;
	  if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
	  break;
	}

	f = atand(sind(asind(tanfi * tant) / 3) / tant);
	hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	if (i > 1 && fabs(swe_difdeg2n(hsp->cusp[ih], cuspsv)) < VERY_SMALL_PLAC_ITER)
	  break;
	cuspsv = hsp->cusp[ih];
      }
      if (i >= niter_max) {
	retc = ERR;
	strcpy(hsp->serr, "very close to polar circle, switched to Porphyry");
	goto porphyry;
      }
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = AscDash(rectasc, f, sine, cose);
#ifdef DEBUG_PLAC_ITER
  fprintf(stderr, "h=%d, niter=%d\n", ih, i);
#endif
    }

    rectasc = swe_degnorm(60 + th);
    tant = tand(asind(sine*sind(Asc1(rectasc,  fh2, sine, cose))));
    ih = 12;
    if (fabs(tant) < VERY_SMALL) {
      hsp->cusp[ih] = rectasc;
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
    } else {
      f = atand(sind(asind(tanfi * tant) / 1.5) / tant);

      hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
      cuspsv = 0;
      for (i = 1; i <= niter_max; i++) {
	tant = tand(asind(sine * sind(hsp->cusp[ih])));
	if (fabs(tant) < VERY_SMALL) {
	  hsp->cusp[ih] = rectasc;
	  if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
	  break;
	}
	f = atand(sind(asind(tanfi * tant) / 1.5) / tant);

	hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	if (i > 1 && fabs(swe_difdeg2n(hsp->cusp[ih], cuspsv)) < VERY_SMALL_PLAC_ITER)
	  break;
	cuspsv = hsp->cusp[ih];
      }
      if (i >= niter_max) {
	retc = ERR;
	strcpy(hsp->serr, "very close to polar circle, switched to Porphyry");
	goto porphyry;
      }
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = AscDash(rectasc, f, sine, cose);
#ifdef DEBUG_PLAC_ITER
  fprintf(stderr, "h=%d, niter=%d\n", ih, i);
#endif
    }

    rectasc = swe_degnorm(120 + th);
    tant = tand(asind(sine * sind(Asc1(rectasc, fh2, sine, cose))));
    ih = 2;
    if (fabs(tant) < VERY_SMALL) {
      hsp->cusp[ih] = rectasc;
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
    } else {
      f = atand(sind(asind(tanfi * tant) / 1.5) / tant);

      hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
      cuspsv = 0;
      for (i = 1; i <= niter_max; i++) {
	tant = tand(asind(sine * sind(hsp->cusp[ih])));
	if (fabs(tant) < VERY_SMALL) {
	  hsp->cusp[ih] = rectasc;
	  if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
	  break;
	}
	f = atand(sind(asind(tanfi * tant) / 1.5) / tant);

	hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	if (i > 1 && fabs(swe_difdeg2n(hsp->cusp[ih], cuspsv)) < VERY_SMALL_PLAC_ITER)
	  break;
	cuspsv = hsp->cusp[ih];
      }
      if (i >= niter_max) {
	retc = ERR;
	strcpy(hsp->serr, "very close to polar circle, switched to Porphyry");
	goto porphyry;
      }
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = AscDash(rectasc, f, sine, cose);
#ifdef DEBUG_PLAC_ITER
  fprintf(stderr, "h=%d, niter=%d\n", ih, i);
#endif
    }

    rectasc = swe_degnorm(150 + th);
    tant = tand(asind(sine * sind(Asc1(rectasc, fh1, sine, cose))));
    ih = 3;
    if (fabs(tant) < VERY_SMALL) {
      hsp->cusp[ih] = rectasc;
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
    } else {
      f = atand(sind(asind(tanfi * tant) / 3) / tant);

      hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
      cuspsv = 0;
      for (i = 1; i <= niter_max; i++) {
	tant = tand(asind(sine * sind(hsp->cusp[ih])));
	if (fabs(tant) < VERY_SMALL) {
	  hsp->cusp[ih] = rectasc;
	  if (hsp->do_hspeed) hsp->cusp_speed[ih] = hsp->armc_speed;
	  break;
	}
	f = atand(sind(asind(tanfi * tant) / 3) / tant);

	hsp->cusp[ih] = Asc1(rectasc, f, sine, cose);
	if (i > 1 && fabs(swe_difdeg2n(hsp->cusp[ih], cuspsv)) < VERY_SMALL_PLAC_ITER)
	  break;
	cuspsv = hsp->cusp[ih];
      }
      if (i >= niter_max) {
	retc = ERR;
	strcpy(hsp->serr, "very close to polar circle, switched to Porphyry");
	goto porphyry;
      }
      if (hsp->do_hspeed) hsp->cusp_speed[ih] = AscDash(rectasc, f, sine, cose);
#ifdef DEBUG_PLAC_ITER
  fprintf(stderr, "h=%d, niter=%d\n", ih, i);
#endif
    }
    break;
  }
  if (hsy != 'G' && hsy != 'Y' && toupper(hsy) != 'I' ) {
    hsp->cusp[4] = swe_degnorm(hsp->cusp[10] + 180);
    hsp->cusp[5] = swe_degnorm(hsp->cusp[11] + 180);
    hsp->cusp[6] = swe_degnorm(hsp->cusp[12] + 180);
    hsp->cusp[7] = swe_degnorm(hsp->cusp[1] + 180);
    hsp->cusp[8] = swe_degnorm(hsp->cusp[2] + 180);
    hsp->cusp[9] = swe_degnorm(hsp->cusp[3] + 180);
    if (hsp->do_hspeed && ! hsp->do_interpol) {
      hsp->cusp_speed[4] = hsp->cusp_speed[10];
      hsp->cusp_speed[5] = hsp->cusp_speed[11];
      hsp->cusp_speed[6] = hsp->cusp_speed[12];
      hsp->cusp_speed[7] = hsp->cusp_speed[1];
      hsp->cusp_speed[8] = hsp->cusp_speed[2];
      hsp->cusp_speed[9] = hsp->cusp_speed[3];
    }
  }

  if (fi >= 0)
    f = 90 - fi;
  else
    f = -90 - fi;
  hsp->vertex = Asc1(th - 90, f, sine, cose);
  if (hsp->do_speed) hsp->vertex_speed = AscDash(th - 90, f, sine, cose);

  if (fabs(fi) <= ekl) {
	vemc = swe_difdeg2n(hsp->vertex, hsp->mc);
    if (vemc > 0)
      hsp->vertex = swe_degnorm(hsp->vertex + 180);
  }

  th2 = swe_degnorm(th + 90);
  if (fabs(th2 - 90) > VERY_SMALL
    && fabs(th2 - 270) > VERY_SMALL) {
    tant = tand(th2);
    hsp->equasc = atand(tant / cose);
    if (th2 > 90 && th2 <= 270)
      hsp->equasc = swe_degnorm(hsp->equasc + 180);
  } else {
    if (fabs(th2 - 90) <= VERY_SMALL)
      hsp->equasc = 90;
    else
      hsp->equasc = 270;
  }
  hsp->equasc = swe_degnorm(hsp->equasc);
  if (hsp->do_speed) hsp->equasc_speed = AscDash(th + 90, 0, sine, cose);

  hsp->coasc1 = swe_degnorm(Asc1(th - 90, fi, sine, cose) + 180);
  if (hsp->do_speed) hsp->coasc1_speed = AscDash(th - 90, fi, sine, cose);

  if (fi >= 0) {
    hsp->coasc2 = Asc1(th + 90, 90 - fi, sine, cose);
    if (hsp->do_speed) hsp->coasc2_speed = AscDash(th + 90, 90 - fi, sine, cose);
  } else {
    hsp->coasc2 = Asc1(th + 90, -90 - fi, sine, cose);
    if (hsp->do_speed) hsp->coasc2_speed = AscDash(th + 90, -90 - fi, sine, cose);
  }

  hsp->polasc = Asc1(th - 90, fi, sine, cose);
  if (hsp->do_speed) hsp->polasc_speed = AscDash(th - 90, fi, sine, cose);
  return retc;
}

static double Asc1(double x1, double f, double sine, double cose)
{
  int n;
  double ass;
  x1 = swe_degnorm(x1);
  n  = (int) ((x1 / 90) + 1);
  if (fabs(90 - f) < VERY_SMALL) {
    return 180;
  }
  if (fabs(90 + f) < VERY_SMALL) {
    return 0;
  }
  if (n == 1)
    ass = ( Asc2(x1, f, sine, cose));
  else if (n == 2)
    ass = (180 - Asc2(180 - x1, - f, sine, cose));
  else if (n == 3)
    ass = (180 + Asc2(x1 - 180, - f, sine, cose));
  else
    ass = (360 - Asc2(360- x1,  f, sine, cose));
  ass = swe_degnorm(ass);
  if (fabs(ass - 90) < VERY_SMALL)
	ass = 90;
  if (fabs(ass - 180) < VERY_SMALL)
    ass = 180;
  if (fabs(ass - 270) < VERY_SMALL)
    ass = 270;
  if (fabs(ass - 360) < VERY_SMALL)
    ass = 0;
  return ass;
}

static double Asc2(double x, double f, double sine, double cose)
{
  double ass, sinx;

  ass = - tand(f) * sine + cose * cosd(x);
  if (fabs(ass) < VERY_SMALL)
    ass = 0;
  sinx = sind(x);
  if (fabs(sinx) < VERY_SMALL)
    sinx = 0;
  if (sinx == 0) {
    if (ass < 0)
      ass = -VERY_SMALL;
    else
      ass = VERY_SMALL;
  } else if (ass == 0) {
    if (sinx < 0)
      ass = -90;
    else
      ass = 90;
  } else {

    ass = atand(sinx / ass);
  }
  if (ass < 0)
    ass = 180 + ass;
  return (ass);
}

static double AscDash(double x, double f, double sine, double cose)
{
  double cosx = cosd(x);
  double sinx = sind(x);
  double sinx2 = sinx * sinx;
  double c = cose * cosx - tand(f) * sine;
  double d = sinx2 + c * c;
  double dudt;
  if (d > VERY_SMALL) {
      dudt = (cosx * c + cose * sinx2) / d;
  } else {
      dudt = 0.0;
  }
  return dudt * ARMCS;
}

static double armc_to_mc(double armc, double eps)
{
  double cose = cosd(eps);
  double mc, tant;
  if (fabs(armc - 90) > VERY_SMALL
	  && fabs(armc - 270) > VERY_SMALL) {
    tant = tand(armc);
    mc = swe_degnorm(atand(tant / cose));
    if (armc > 90 && armc <= 270)
    mc = swe_degnorm(mc + 180);
  } else {
    if (fabs(armc - 90) <= VERY_SMALL)
      mc = 90;
    else
      mc = 270;
  }
  return mc;
}

static double fix_asc_polar(double asc, double armc, double eps, double geolat)
{
  double demc = atand(sind(armc) * tand(eps));
  if (geolat >= 0 && 90 - geolat + demc < 0)
    asc = swe_degnorm(asc + 180);
  if (geolat < 0 && -90 - geolat + demc > 0)
    asc = swe_degnorm(asc + 180);
  return asc;
}

double CALL_CONV swe_house_pos(
	double armc, double geolat, double eps, int hsys, double *xpin, char *serr)
{
  double xp[6], xeq[6], ra, de, mdd, mdn, sad, san;
  double hpos, sinad, ad, a, admc, adp, samc, asc, mc, acmc, tant;

  double fh, ra0, tanfi, sinfi, fac, dfac, tanx;
  double x[3], xasc[3], xs1, xs2, raep, raaz, oblaz, xtemp;
  double hcusp[37], ascmc[10];
  double sine = sind(eps);
  double cose = cosd(eps);
  double c1, c2 = 0.0, d, hsize;
  int i, j, nloop;
  double dsun = 0, darmc, harmc, y, sinpsi, sa;
  AS_BOOL is_western_half = FALSE;
  hsys = toupper(hsys);
  if (1) {

    ascmc[9] = 99;

    if (swe_houses_armc_ex2(armc, geolat, eps, hsys, hcusp, ascmc, NULL, NULL, serr) == ERR) {
      if (serr != NULL)
	sprintf(serr, "swe_house_pos(): failed for system %c", hsys);
    } else {
      hpos = 0;
      for (i = 1; i <= 12; i++) {
	if (fabs(swe_difdeg2n(xpin[0], hcusp[i])) < MILLIARCSEC && xpin[1] == 0) {
	  hpos = (double) i;
	}
      }
      for (i = 1; i <= 12; i += 3) {
	if (fabs(swe_difdeg2n(xpin[0], hcusp[i])) < MILLIARCSEC && xpin[1] == 0) {
	  xp[0] = (double) i;
	}
      }
      if (hpos > 0)
	return hpos;

      if (hsys == 'I')
	dsun = ascmc[9];

      if (hsys == 'Y') {
	xeq[0] = ascmc[0];
	xeq[1] = 0;
	xeq[2] = 1;
	swe_cotrans(xeq, xeq, -eps);
	dsun = xeq[1];
      }
    }
  }
  AS_BOOL is_above_hor = FALSE;
  AS_BOOL is_invalid = FALSE;
  AS_BOOL is_circumpolar = FALSE;
  if (serr != NULL)
    *serr = '\0';
  xeq[0] = xpin[0];
  xeq[1] = xpin[1];
  xeq[2] = 1;
  swe_cotrans(xeq, xeq, -eps);
  ra = xeq[0];
  de = xeq[1];
  mdd = swe_degnorm(ra - armc);
  mdn = swe_degnorm(mdd + 180);
  if (mdd >= 180)
    mdd -= 360;
  if (mdn >= 180)
    mdn -= 360;

  switch(hsys) {
    case 'N':
      xp[0] = xpin[0];
      hpos = xp[0] / 30.0 + 1;
      break;
    case 'A':
    case 'E':
    case 'D':
    case 'V':
    case 'W':
      asc = Asc1(swe_degnorm(armc + 90), geolat, sine, cose);
      mc = armc_to_mc(armc, eps);
      asc = fix_asc_polar(asc, armc, eps, geolat);
      xp[0] = swe_degnorm(xpin[0] - asc);
      if (hsys == 'V')
	xp[0] = swe_degnorm(xp[0] + 15);
      if (hsys == 'W')
	xp[0] = swe_degnorm(xp[0] + fmod(asc, 30));
      if (hsys == 'D')
	xp[0] = swe_degnorm(xpin[0] - mc - 90);

      xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      hpos = xp[0] / 30.0 + 1;
      break;
    case 'O':
    case 'B':
    case 'S':
      asc = Asc1(swe_degnorm(armc + 90), geolat, sine, cose);

      mc = armc_to_mc(armc, eps);

      asc = fix_asc_polar(asc, armc, eps, geolat);
      if (hsys ==  'O' || hsys == 'S') {
	xp[0] = swe_degnorm(xpin[0] - asc);

	xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
	if (xp[0] < 180)
	  hpos = 1;
	else {
	  hpos = 7;
	  xp[0] -= 180;
	}
	acmc = swe_difdeg2n(asc, mc);
	if (xp[0] < 180 - acmc)
	  hpos += xp[0] * 3 / (180 - acmc);
	else
	  hpos += 3 + (xp[0] - 180 + acmc) * 3 / acmc;
        if (hsys == 'S') {
	  hpos += 0.5;
	  if (hpos > 12) hpos = 1;
	}
      } else {
	double dek, r, sna, sda;
	dek = asind(sind(asc) * sine);

	tanfi = tand(geolat);
	r = -tanfi * tand(dek);

	sda = acos(r) * RADTODEG;
	sna = 180 - sda;
	if (mdd > 0) {
	  if (mdd < sda)
	    hpos = mdd * 90 / sda;
	  else
	    hpos = 90 + (mdd - sda) * 90 / sna;
	} else {
	  if (mdd > -sna)
	    hpos = 360 + mdd * 90 / sna;
	  else
	    hpos = 270 + (mdd + sna) * 90 / sda;
	}
	hpos = swe_degnorm(hpos - 90) / 30.0 + 1.0;
	if (hpos >= 13.0) hpos -= 12;
      }
      break;
    case 'X':
      hpos = swe_degnorm(mdd - 90) / 30.0 + 1.0;
      break;
    case 'F':
      x[0] = Asc1(swe_degnorm(armc + 90), geolat, sine, cose);
      x[0] = fix_asc_polar(x[0], armc, eps, geolat);
      x[1] = 0;
      swe_cotrans(x, x, -eps);
      hpos = swe_degnorm(ra - x[0]) / 30.0 + 1;
      break;
    case 'M': {
      double a = xpin[0];
      if (fabs(a - 90) > VERY_SMALL
        && fabs(a - 270) > VERY_SMALL) {
        tant = tand(a);
	hpos = atand(tant / cose);
        if (a > 90 && a <= 270)
          hpos = swe_degnorm(hpos + 180);
      } else {
	if (fabs(a - 90) <= VERY_SMALL)
          hpos = 90;
        else
          hpos = 270;
      }
      hpos = swe_degnorm(hpos - armc - 90);
      hpos = hpos / 30.0 + 1;
    }
      break;

    case 'K':

      is_invalid = FALSE;
      is_circumpolar = FALSE;

      if (90 - geolat < de || -90 - geolat > de) {
        adp = 90;
	is_circumpolar = TRUE;
      }

      else if (geolat - 90 > de || geolat + 90 < de) {
        adp = -90;
	is_circumpolar = TRUE;
      }

      else {
	adp = asind(tand(geolat) * tand(de));
      }
      admc = tand(eps) * tand(geolat) * sind(armc);

      if (fabs(admc) > 1) {
	if (admc > 1)
	  admc = 1;
	else
	  admc = -1;
	is_circumpolar = TRUE;
      }
      admc = asind(admc);
      samc = 90 + admc;
      if (samc == 0)
        is_invalid = TRUE;
      if (fabs(samc) > 0) {
	if (mdd >= 0) {
	  dfac = (mdd - adp + admc) / samc;
	  xp[0] = swe_degnorm((dfac - 1) * 90);
	  xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);

	  if (dfac > 2 || dfac < 0)
	    is_invalid = TRUE;
	} else {
	  dfac = (mdd + 180 + adp + admc) / samc;
	  xp[0] = swe_degnorm((dfac + 1) * 90);
	  xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);

	  if (dfac > 2 || dfac < 0)
	    is_invalid = TRUE;
	}
      }
      if (is_invalid) {
        xp[0] = 0;
	hpos = 0;
	if (serr != NULL)
          strcpy(serr, "Koch house position failed in circumpolar area");
	break;
      }
      if (is_circumpolar) {
	if (serr != NULL)
          strcpy(serr, "Koch house position, doubtful result in circumpolar area");
      }

      hpos = xp[0] / 30.0 + 1;
      break;
    case 'C':
      xeq[0] = swe_degnorm(mdd - 90);

      swe_cotrans(xeq, xp, -geolat);

      xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      hpos = xp[0] / 30.0 + 1;
      break;
    case 'J':
      sinfi = sind(geolat);
      if (fabs(geolat) < VERY_SMALL) {
	xs2 = 1 / 3.0;
	xs1 = 2 / 3.0;
      } else {
	xs2 = sind(geolat / 3) / sinfi;
	xs1 = sind(2 * geolat / 3) / sinfi;
      }
      xs2 = asind(xs2);
      xs1 = asind(xs1);

      hcusp[1] = 0;
      hcusp[2] = xs2;
      hcusp[3] = xs1;
      hcusp[4] = 90;
      hcusp[5] = 180 - xs1;
      hcusp[6] = 180 - xs2;
      hcusp[7] = 180;
      hcusp[8] = 180 + xs2;
      hcusp[9] = 180 + xs1;
      hcusp[10] = 270;
      hcusp[11] = 360 - xs1;
      hcusp[12] = 360 - xs2;
      xeq[0] = swe_degnorm(mdd - 90);
      swe_cotrans(xeq, xp, -geolat);
      a = xp[0];
      if (swe_difdeg2n(hcusp[6], hcusp[1]) > 0) {
	d = swe_degnorm(a - hcusp[1]);
	for (i = 1; i <= 12; i++) {
	  j = i + 1;
	  if (j > 12)
	    c2 = 360;
	  else
	    c2 = swe_degnorm(hcusp[j] - hcusp[1]);
	  if (d < c2) break;
	}
	c1 = swe_degnorm(hcusp[i] - hcusp[1]);
      } else {
	d = swe_degnorm(hcusp[1] - a);
	for (i = 1; i <= 12; i++) {
	  j = i + 1;
	  if (j > 12)
	    c2 = 360;
	  else
	    c2 = swe_degnorm(hcusp[1] - hcusp[j]);
	  if (d < c2) break;
	}
	c1 = swe_degnorm(hcusp[1] - hcusp[i]);
      }
      hsize = c2 - c1;
      if (hsize == 0) {
	hpos = i;
      } else {
	hpos = i + (d - c1) / hsize;
      }
      break;
    case 'U':
      if (fabs(geolat) < VERY_SMALL) {
        geolat = (geolat >= 0) ? VERY_SMALL : -VERY_SMALL;
      }

      asc = Asc1(swe_degnorm(armc + 90), geolat, sine, cose);

      asc = fix_asc_polar(asc, armc, eps, geolat);

      x[0] = asc; x[1] = 0.0; x[2] = 1.0;
      swe_cotrans(x, x, -eps);
      raep = swe_degnorm(armc + 90);
      x[0] = swe_degnorm(raep - x[0]);
      swe_cotrans(x, x, -(90-geolat));
      tanx = tand(x[0]);
      if (geolat == 0) {
        xtemp = (tanx >= 0) ? 90 : -90;
      } else {
	xtemp = atand(tanx/cosd((90-geolat)));
      }
      if (x[0] > 90 && x[0] <= 270)
	xtemp = swe_degnorm(xtemp + 180);
      x[0] = swe_degnorm(xtemp);
      raaz = swe_degnorm(raep - x[0]);

      x[0] = raaz; x[1] = 0.0;
      x[0] = swe_degnorm(raep - x[0]);
      swe_cotrans(x, x, -(90-geolat));
      x[1] = x[1] + 90;
      swe_cotrans(x, x, 90-geolat);
      oblaz = x[1];

      xasc[0] = asc; xasc[1] = 0.0; xasc[2] = 1.0;
      swe_cotrans(xasc, xasc, -eps);
      xasc[0] = swe_degnorm(xasc[0] - raaz);
      xtemp = atand(tand(xasc[0])/cosd(oblaz));
      if (xasc[0] > 90 && xasc[0] <= 270)
          xtemp = swe_degnorm(xtemp + 180);
      xasc[0] = swe_degnorm(xtemp);

      xp[0] = swe_degnorm(xeq[0] - raaz);
      xtemp = atand(tand(xp[0])/cosd(oblaz));
      if (xp[0] > 90 && xp[0] <= 270)
	xtemp = swe_degnorm(xtemp + 180);
      xp[0] = swe_degnorm(xtemp);
      xp[0] = swe_degnorm(xp[0]-xasc[0]);

      x[0] = xeq[0];
      x[1] = xeq[1];
      swe_cotrans(x, x, oblaz);
      xp[1] = xeq[1] - x[1];

      xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      hpos = xp[0] / 30.0 + 1;
      break;
    case 'H':
      xeq[0] = swe_degnorm(mdd - 90);
      swe_cotrans(xeq, xp, 90 - geolat);

      xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      hpos = xp[0] / 30.0 + 1;
      break;
    case 'R':
      if (fabs(mdd) < VERY_SMALL)
	xp[0] = 270;
      else if (180 - fabs(mdd) < VERY_SMALL)
        xp[0] = 90;
      else {
        if (90 - fabs(geolat) < VERY_SMALL) {
          if (geolat > 0)
	    geolat = 90 - VERY_SMALL;
          else
	    geolat = -90 + VERY_SMALL;
        }
        if (90 - fabs(de) < VERY_SMALL) {
          if (de > 0)
            de = 90 - VERY_SMALL;
          else
	    de = -90 + VERY_SMALL;
        }
        a = tand(geolat) * tand(de) + cosd(mdd);
        xp[0] = swe_degnorm(atand(-a / sind(mdd)));
        if (mdd < 0)
          xp[0] += 180;
        xp[0] = swe_degnorm(xp[0]);

        xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      }
      hpos = xp[0] / 30.0 + 1;
      break;

    case 'I': case 'i':
    case 'Y':
      if (geolat > 90 - MILLIARCSEC)
        geolat = 90 - MILLIARCSEC;
      if (geolat < -90 + MILLIARCSEC)
        geolat = -90 + MILLIARCSEC;

      if (90 - fabs(de) < VERY_SMALL) {
	if (de > 0)
	  de = 90 - VERY_SMALL;
	else
	  de = -90 + VERY_SMALL;
      }
      a = tand(geolat) * tand(de) + cosd(mdd);
      xp[0] = swe_degnorm(atand(-a / sind(mdd)));
      if (mdd < 0)
	xp[0] += 180;
      xp[0] = swe_degnorm(xp[0]);

      sinad = tand(de) * tand(geolat);
      a = sinad + cosd(mdd);
      if (a >= 0)
	is_above_hor = TRUE;

      harmc = 90 - geolat;
      if (geolat < 0)
	harmc = 90 + geolat;

      darmc = swe_degnorm(xp[0] - 270);
      if (darmc > 180) {
	is_western_half = TRUE;
	darmc = (360 - darmc);
      }

      sinad = tand(dsun) * tand(geolat);
      if (sinad >= 1)
	ad = 90;

      else if (sinad <= -1)
	ad = -90;

      else
	ad = asind(sinad);
      sad = 90 + ad;
      san = 90 - ad;

      if (sad == 0 && is_above_hor) {
	xp[0] = 270;

      } else if (san == 0 && !is_above_hor) {
	xp[0] = 90;

      } else {
	sa = sad;
	if (!is_above_hor) {
	  dsun = -dsun;
	  sa = san;
	  darmc = 180 - darmc;
	  is_western_half = !is_western_half;
	}

	a = acosd(cosd(harmc) * cosd(darmc));
	if (a < VERY_SMALL)
	  a = VERY_SMALL;

	sinpsi = sind(harmc) / sind(a);
	if (sinpsi > 1) sinpsi = 1;
	if (sinpsi < -1) sinpsi = -1;

	y = sind(dsun) / sinpsi;
	if (y > 1)
	  y = 90 - VERY_SMALL;
	else if (y < -1)
	  y = - (90 - VERY_SMALL);
	else
	  y = asind(y);
	d = acosd(cosd(y) / cosd(dsun));
	if (dsun < 0) d = -d;
	if (geolat < 0) d = -d;
	darmc += d;
	if (is_western_half)
	  xp[0] = 270 - (darmc / sa) * 90;
	else
	  xp[0] = 270 + (darmc / sa) * 90;
	if (!is_above_hor)
	  xp[0] = swe_degnorm(xp[0] + 180);
      }

      xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      hpos = xp[0] / 30.0 + 1;
      break;
    case 'T':
      fh = geolat;
      if (fh > 89.999)
	fh = 89.999;
      if (fh < -89.999)
	fh = -89.999;
      mdd = swe_degnorm(mdd);
      if (de > 90 - VERY_SMALL)
	de = 90 - VERY_SMALL;
      if (de < -90 + VERY_SMALL)
	de = -90 + VERY_SMALL;
      sinad = tand(de) * tand(fh);
      if (sinad > 1.0) sinad = 1.0;
      if (sinad < -1.0) sinad = -1.0;
      a = sinad + cosd(mdd);
      if (a >= 0)
        is_above_hor = TRUE;

      if (!is_above_hor) {
	ra = swe_degnorm(ra + 180);
	de = -de;
	mdd = swe_degnorm(mdd + 180);
      }

      if (mdd > 180) {
	ra = swe_degnorm(armc - mdd);
      }

      tanfi = tand(fh);
      ra0 = swe_degnorm(armc + 90);
      xp[1] = 1;
      xeq[1] = de;
      fac = 2;
      nloop = 0;
      while (fabs(xp[1]) > 0.000001 && nloop < 1000) {
	if (xp[1] > 0) {
	  fh = atand(tand(fh) - tanfi / fac);
	  ra0 -= 90 / fac;
	} else {
	  fh = atand(tand(fh) + tanfi / fac);
	  ra0 += 90 / fac;
	}
	xeq[0] = swe_degnorm(ra - ra0);
        swe_cotrans(xeq, xp, 90 - fh);
	fac *= 2;
	nloop++;
      }
      hpos = swe_degnorm(ra0 - armc);

      if (mdd > 180)
	    hpos = swe_degnorm(-hpos);

      if (!is_above_hor)
	hpos = swe_degnorm(hpos + 180);
      hpos = swe_degnorm(hpos - 90) / 30 + 1;
      break;
    case 'P':
    case 'G':

      if (90 - fabs(de) <= fabs(geolat)) {
        if (de * geolat < 0)
          xp[0] = swe_degnorm(90 + mdn / 2);
        else
          xp[0] = swe_degnorm(270 + mdd / 2);
	if (serr != NULL)
          strcpy(serr, "Otto Ludwig procedure within circumpolar regions.");
      } else {
        sinad = tand(de) * tand(geolat);
        ad = asind(sinad);
        a = sinad + cosd(mdd);
        if (a >= 0)
          is_above_hor = TRUE;
        sad = 90 + ad;
        san = 90 - ad;
        if (is_above_hor)
          xp[0] =  (mdd / sad + 3) * 90;
        else
          xp[0] = (mdn / san + 1) * 90;

        xp[0] = swe_degnorm(xp[0] + MILLIARCSEC);
      }
      if (hsys == 'G') {
        xp[0] = 360 - xp[0];
        hpos = xp[0] / 10.0 + 1;
      } else {
        hpos = xp[0] / 30.0 + 1;
      }
    break;
  default:
    hpos = 0;
    if (swe_houses_armc_ex2(armc, geolat, eps, hsys, hcusp, ascmc, NULL, NULL, serr) == ERR) {
      if (serr != NULL)
	sprintf(serr, "swe_house_pos(): failed for system %c", hsys);
      break;
    }
    if (swe_difdeg2n(hcusp[6], hcusp[1]) > 0) {
      d = swe_degnorm(xpin[0] - hcusp[1]);
      for (i = 1; i <= 12; i++) {
	j = i + 1;
	if (j > 12)
	  c2 = 360;
	else
	  c2 = swe_degnorm(hcusp[j] - hcusp[1]);
	if (d < c2) break;
      }
      c1 = swe_degnorm(hcusp[i] - hcusp[1]);
    } else {
      d = swe_degnorm(hcusp[1] - xpin[0]);
      for (i = 1; i <= 12; i++) {
	j = i + 1;
	if (j > 12)
	  c2 = 360;
	else
	  c2 = swe_degnorm(hcusp[1] - hcusp[j]);
	if (d < c2) break;
      }
      c1 = swe_degnorm(hcusp[1] - hcusp[i]);
    }
    hsize = c2 - c1;
    if (hsize == 0) {
      hpos = i;
    } else {
      hpos = i + (d - c1) / hsize;
    }
    if (serr != NULL)
      sprintf(serr, "swe_house_pos(): using simplified algorithm for system %c\n", hsys);
    break;
  }
  return hpos;
}

static int sunshine_init(double lat, double dec, double xh[])
{
  double ad, nsa, dsa, arg;

  arg = tand(dec) * tand(lat);
  if (arg >= 1) {
    ad = 90 - VERY_SMALL;
  } else if (arg <= -1) {
    ad = -90 + VERY_SMALL;
  } else {
    ad = asind(arg);
  }
  nsa = 90 - ad;
  dsa = 90 + ad;
  xh[2] = -2 * nsa / 3;
  xh[3] = -1 * nsa / 3;
  xh[5] = 1 * nsa / 3;
  xh[6] = 2 * nsa / 3;
  xh[8] = -2 * dsa / 3;
  xh[9] = -1 * dsa / 3;
  xh[11] = 1 * dsa / 3;
  xh[12] = 2 * dsa / 3;
  if (fabs(arg) >= 1)
    return ERR;
  return OK;
}

static int sunshine_solution_makransky(double ramc, double lat, double ecl, struct houses *hsp)
{
  double xh[13];
  double md;
  double zd;
  double pole, q, w, a, b, c, f, cu, r = 0, rah;
  double sinlat, coslat, tanlat, tandec, sinecl;
  double dec = hsp->sundec;
  sinlat = sind(lat);
  coslat = cosd(lat);
  tanlat = tand(lat);
  tandec = tand(dec);
  sinecl = sind(ecl);
  int ih;

  if (sunshine_init(lat, dec, xh) == ERR)
    return ERR;
  for (ih = 1; ih <= 12; ih++) {
    double z = 0;
    if ((ih - 1) % 3 == 0) continue;
    md = fabs(xh[ih]);
    if (ih <= 6)
      rah = swe_degnorm(ramc + 180 + xh[ih]);
    else
      rah = swe_degnorm(ramc + xh[ih]);
    if (lat < 0) {
      rah = swe_degnorm(180 + rah);
    }

    if (md == 90) {

      zd = 90.0 - atand(sinlat * tandec);
    } else {
      if (md < 90) {

	a = atand(coslat * tand(md));
      } else {

	a = atand(tand(md - 90) / coslat);
      }

      b = atand(tanlat * cosd(md));

      if (ih <= 6)
	c = b + dec;
      else
	c = b - dec;

      f = atand(sinlat * sind(md) * tand(c));

      zd = a + f;
    }
    pole = asind(sind(zd) * sinlat);
    q = asind(tandec * tand(pole));
    if (ih <= 3 || ih >= 11)
      w = swe_degnorm(rah - q);
    else
      w = swe_degnorm(rah + q);
    if (w == 90) {
      r = atand(sind(ecl) * tand(pole));
      if (ih <= 3 || ih >= 11)
        cu = 90 + r;
      else
        cu = 90 - r;
    } else if (w == 270) {
      r = atand(sinecl * tand(pole));
      if (ih <= 3 || ih >= 11)
        cu = 270 - r;
      else
        cu = 270 + r;
    } else {
      double m;
      m = atand(fabs(tand(pole) / cosd(w)));
      if (ih <= 3 || ih >= 11) {
        if (w > 90 && w < 270)
	  z = m - ecl;
	else
	  z = m + ecl;
      } else {
        if (w > 90 && w < 270)
	  z = m + ecl;
	else
	  z = m - ecl;
      }
      if (z == 90) {
        if (w < 180)
	  cu = 90;
	else
	  cu = 270;
      } else {

        r = atand(fabs(cosd(m) * tand(w) / cosd(z)));
	if (w < 90)
	  cu = r;
	else if (w > 90 && w < 180)
	  cu = 180 - r;
	else if (w > 180 && w < 270)
	  cu = 180 + r;
	else
	  cu = 360 - r;
      }
      if (z > 90) {

	if (w < 90)
	  cu = 180 - r;
	else if (w > 90 && w < 180)
	  cu = + r;
	else if (w > 180 && w < 270)
	  cu = 360 - r;
	else
	  cu = 180 + r;
      }
      if (lat < 0)
        cu = swe_degnorm(cu + 180);
    }
    hsp->cusp[ih] = cu;
  }
  return OK;
}

static int sunshine_solution_treindl(double ramc, double lat, double ecl, struct houses *hsp)
{
  double xh[13];
  double mcdec, sinlat, coslat, cosdec, tandec, sinecl, cosecl;
  double xhs, pole, a, cosa, alph, alpha2, c, cosc, b, sinzd, zd, rax, hc;
  int ih, retval = OK;
  AS_BOOL mc_under_horizon;
  double dec = hsp->sundec;

  sinlat = sind(lat);
  coslat = cosd(lat);
  cosdec = cosd(dec);
  tandec = tand(dec);
  sinecl = sind(ecl);
  cosecl = cosd(ecl);
  sunshine_init(lat, dec, xh);

  mcdec = atand(sind(ramc) * tand(ecl));
  mc_under_horizon = fabs(lat - mcdec) > 90;
  if (mc_under_horizon && SUNSHINE_KEEP_MC_SOUTH) {

    for (ih = 2; ih <= 12; ih++) {
      xh[ih] = -xh[ih];
    }
  }

  for (ih = 1; ih <= 12; ih++) {
    if ((ih - 1) % 3 == 0) continue;
    xhs = 2 * asind(cosdec * sind(xh[ih] / 2));

    cosa = tandec * tand(xhs / 2);
    alph = acosd(cosa);

    if (ih > 7) {

      alpha2 = 180 - alph;
      b = 90 - lat + dec;
    } else {
      alpha2 = alph;
      b = 90 - lat - dec;
    }

    cosc = cosd(xhs) * cosd(b) + sind(xhs) * sind(b) * cosd(alpha2);
    c = acosd(cosc);

    if (c < 1e-6) {
      sprintf(hsp->serr, "Sunshine house %d c=%le very small", ih, c);
      retval = ERR;
    }
    sinzd = sind(xhs) * sind(alpha2) / sind(c);
    zd = asind(sinzd);

    rax = atand(coslat * tand(zd));

    pole = asind(sinzd * sinlat);
    if (ih <= 6) {
      pole = -pole;
      a = swe_degnorm(rax + ramc + 180);
    } else {
      a = swe_degnorm(ramc + rax);
    }

    hc = Asc1(a, pole, sinecl, cosecl);
    hsp->cusp[ih] = hc;
  }
  if (mc_under_horizon && ! SUNSHINE_KEEP_MC_SOUTH) {
    for (ih = 2; ih <= 12; ih++) {
      if ((ih - 1) % 3 == 0) continue;
      hsp->cusp[ih] = swe_degnorm(hsp->cusp[ih] + 180);
    }
  }
  return retval;
}
