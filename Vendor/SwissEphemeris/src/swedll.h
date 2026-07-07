#ifdef __cplusplus
extern "C" {
#endif
#ifndef _SWEDLL_H
#define _SWEDLL_H

#ifndef _SWEPHEXP_INCLUDED
#include "swephexp.h"
#endif

# ifdef __cplusplus
#define DllImport extern "C" __declspec( dllimport )
# else
#define DllImport  __declspec( dllimport )
# endif

#if defined (PASCAL) || defined(__stdcall)
  #if defined UNDECO_DLL
    #define CALL_CONV_IMP __cdecl
  #else
    #define CALL_CONV_IMP __stdcall
  #endif
#else
  #define CALL_CONV_IMP
#endif

DllImport int32 CALL_CONV_IMP swe_heliacal_ut(double JDNDaysUTStart, double *geopos, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 iflag, double *dret, char *serr);
DllImport int32 CALL_CONV_IMP swe_heliacal_pheno_ut(double JDNDaysUT, double *geopos, double *datm, double *dobs, char *ObjectName, int32 TypeEvent, int32 helflag, double *darr, char *serr);
DllImport int32 CALL_CONV_IMP swe_vis_limit_mag(double tjdut, double *geopos, double *datm, double *dobs, char *ObjectName, int32 helflag, double *dret, char *serr);

DllImport int32 CALL_CONV_IMP swe_heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr);
DllImport int32 CALL_CONV_IMP swe_topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr);

DllImport double CALL_CONV_IMP swe_degnorm(double deg);

DllImport char * CALL_CONV_IMP swe_version(char *);
DllImport char * CALL_CONV_IMP swe_get_library_path(char *);

DllImport int32 CALL_CONV_IMP swe_calc(
        double tjd, int ipl, int32 iflag,
        double *xx,
        char *serr);
DllImport int32 CALL_CONV_IMP  swe_calc_pctr(
        double tjd, int32 ipl, int32 iplctr, int32 iflag,
	double *xxret,
	char *serr);

DllImport int32 CALL_CONV_IMP swe_calc_ut(
        double tjd_ut, int32 ipl, int32 iflag,
        double *xx,
        char *serr);

DllImport double CALL_CONV_IMP swe_solcross(
	double x2cross, double jd_et, int32 flag, char *serr);
DllImport double CALL_CONV_IMP swe_solcross_ut(
	double x2cross, double jd_ut, int32 flag, char *serr);
DllImport double CALL_CONV_IMP swe_mooncross(
	double x2cross, double jd_et, int32 flag, char *serr);
DllImport double CALL_CONV_IMP swe_mooncross_ut(
	double x2cross, double jd_ut, int32 flag, char *serr);
DllImport double CALL_CONV_IMP swe_mooncross_node(
	double jd_et, int32 flag, double *xlon, double *xlat, char *serr);
DllImport double CALL_CONV_IMP swe_mooncross_node_ut(
	double jd_ut, int32 flag, double *xlon, double *xlat, char *serr);
DllImport int32 CALL_CONV_IMP swe_helio_cross(
	int ipl, double x2cross, double jd_et, int32 iflag, int32 dir, double *jd_cross, char *serr);
DllImport int32 CALL_CONV_IMP swe_helio_cross_ut(
	int ipl, double x2cross, double jd_ut, int32 iflag, int32 dir, double *jd_cross, char *serr);

DllImport int32 CALL_CONV_IMP swe_fixstar(
        char *star, double tjd, int32 iflag,
        double *xx,
        char *serr);

DllImport int32 CALL_CONV_IMP swe_fixstar_ut(
        char *star, double tjd_ut, int32 iflag,
        double *xx,
        char *serr);

DllImport int32 CALL_CONV_IMP swe_fixstar_mag(
        char *star, double *xx, char *serr);

DllImport int32 CALL_CONV_IMP swe_fixstar2(
        char *star, double tjd, int32 iflag,
        double *xx,
        char *serr);

DllImport int32 CALL_CONV_IMP swe_fixstar2_ut(
        char *star, double tjd_ut, int32 iflag,
        double *xx,
        char *serr);

DllImport int32 CALL_CONV_IMP swe_fixstar2_mag(
        char *star, double *xx, char *serr);

DllImport double CALL_CONV_IMP swe_sidtime0(double tjd_ut, double ecl, double nut);
DllImport double CALL_CONV_IMP swe_sidtime(double tjd_ut);

DllImport double CALL_CONV_IMP swe_deltat_ex(double tjd, int32 iflag, char *serr);
DllImport double CALL_CONV_IMP swe_deltat(double tjd);

DllImport int  CALL_CONV_IMP swe_houses(
        double tjd_ut, double geolat, double geolon, int hsys,
        double *hcusps, double *ascmc);

DllImport int  CALL_CONV_IMP swe_houses_ex(
        double tjd_ut, int32 iflag, double geolat, double geolon, int hsys,
        double *hcusps, double *ascmc);

DllImport int  CALL_CONV_IMP swe_houses_ex2(
        double tjd_ut, int32 iflag, double geolat, double geolon, int hsys,
        double *hcusps, double *ascmc, double *cusp_speed, double *ascmc_speed, char *serr);

DllImport int  CALL_CONV_IMP swe_houses_armc(
        double armc, double geolat, double eps, int hsys,
        double *hcusps, double *ascmc);

DllImport int  CALL_CONV_IMP swe_houses_armc_ex2(
        double armc, double geolat, double eps, int hsys,
        double *hcusps, double *ascmc, double *cusp_speed, double *ascmc_speed, char *serr);

DllImport double  CALL_CONV_IMP swe_house_pos(
        double armc, double geolon, double eps, int hsys, double *xpin, char *serr);

DllImport const char * CALL_CONV_IMP swe_house_name(int hsys);

DllImport int32  CALL_CONV_IMP swe_gauquelin_sector(
	double t_ut, int32 ipl, char *starname, int32 iflag, int32 imeth, double *geopos, double atpress, double attemp, double *dgsect, char *serr);

DllImport void  CALL_CONV_IMP swe_set_sid_mode(
        int32 sid_mode, double t0, double ayan_t0);

DllImport int32  CALL_CONV_IMP swe_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr);
DllImport int32  CALL_CONV_IMP swe_get_ayanamsa_ex_ut(double tjd_ut, int32 iflag, double *daya, char *serr);

DllImport double  CALL_CONV_IMP swe_get_ayanamsa(double tjd_et);
DllImport double  CALL_CONV_IMP swe_get_ayanamsa_ut(double tjd_ut);

DllImport char * CALL_CONV_IMP swe_get_ayanamsa_name(int32 isidmode);
DllImport char * CALL_CONV_IMP swe_get_current_file_data(int ifno, double *tfstart, double *tfend, int *denum);

DllImport int  CALL_CONV_IMP swe_date_conversion(
        int y , int m , int d ,
        double utime,
        char c,
        double *tjd);

DllImport double  CALL_CONV_IMP swe_julday(
        int year, int mon, int mday,
        double hour,
        int gregflag);

DllImport void  CALL_CONV_IMP swe_revjul(
        double jd, int gregflag,
        int *year, int *mon, int *mday,
        double *hour);

DllImport void  CALL_CONV_IMP swe_utc_time_zone(
        int32 iyear, int32 imonth, int32 iday,
	int32 ihour, int32 imin, double dsec,
	double d_timezone,
	int32 *iyear_out, int32 *imonth_out, int32 *iday_out,
	int32 *ihour_out, int32 *imin_out, double *dsec_out);

DllImport int32  CALL_CONV_IMP swe_utc_to_jd(
        int32 iyear, int32 imonth, int32 iday,
	int32 ihour, int32 imin, double dsec,
	int32 gregflag, double *dret, char *serr);

DllImport void  CALL_CONV_IMP swe_jdet_to_utc(
        double tjd_et, int32 gregflag,
	int32 *iyear, int32 *imonth, int32 *iday,
	int32 *ihour, int32 *imin, double *dsec);

DllImport void  CALL_CONV_IMP swe_jdut1_to_utc(
        double tjd_ut, int32 gregflag,
	int32 *iyear, int32 *imonth, int32 *iday,
	int32 *ihour, int32 *imin, double *dsec);

DllImport int  CALL_CONV_IMP swe_time_equ(
        double tjd, double *e, char *serr);
DllImport int  CALL_CONV_IMP swe_lmt_to_lat(double tjd_lmt, double geolon, double *tjd_lat, char *serr);
DllImport int  CALL_CONV_IMP swe_lat_to_lmt(double tjd_lat, double geolon, double *tjd_lmt, char *serr);

DllImport double  CALL_CONV_IMP swe_get_tid_acc(void);
DllImport void  CALL_CONV_IMP swe_set_tid_acc(double tidacc);
DllImport void  CALL_CONV_IMP swe_set_delta_t_userdef(double dt);
DllImport void  CALL_CONV_IMP swe_set_ephe_path(const char *path);
DllImport void  CALL_CONV_IMP swe_set_jpl_file(const char *fname);
DllImport void  CALL_CONV_IMP swe_close(void);
DllImport char * CALL_CONV_IMP swe_get_planet_name(int ipl, char *spname);
DllImport void  CALL_CONV_IMP swe_cotrans(double *xpo, double *xpn, double eps);
DllImport void  CALL_CONV_IMP swe_cotrans_sp(double *xpo, double *xpn, double eps);

DllImport void  CALL_CONV_IMP swe_set_topo(double geolon, double geolat, double height);

DllImport void CALL_CONV_IMP swe_set_astro_models(char *samod, int32 iflag);
DllImport void CALL_CONV_IMP swe_get_astro_models(char *samod, char *sdet, int32 iflag);

DllImport int32  CALL_CONV_IMP swe_sol_eclipse_where(double tjd, int32 ifl, double *geopos, double *attr, char *serr);

DllImport int32  CALL_CONV_IMP swe_lun_occult_where(double tjd, int32 ipl, char *starname, int32 ifl, double *geopos, double *attr, char *serr);

DllImport int32  CALL_CONV_IMP swe_sol_eclipse_how(double tjd, int32 ifl, double *geopos, double *attr, char *serr);

DllImport int32  CALL_CONV_IMP swe_sol_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr);

DllImport int32  CALL_CONV_IMP swe_lun_occult_when_loc(double tjd_start, int32 ipl, char *starname, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr);

DllImport int32  CALL_CONV_IMP swe_sol_eclipse_when_glob(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr);

DllImport int32  CALL_CONV_IMP swe_lun_occult_when_glob(double tjd_start, int32 ipl, char *starname, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr);

DllImport int32  CALL_CONV_IMP swe_lun_eclipse_how(
          double tjd_ut,
          int32 ifl,
	  double *geopos,
          double *attr,
          char *serr);
DllImport int32  CALL_CONV_IMP swe_lun_eclipse_when(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr);
DllImport int32  CALL_CONV_IMP swe_lun_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr);

DllImport int32  CALL_CONV_IMP swe_pheno(double tjd, int32 ipl, int32 iflag, double *attr, char *serr);

DllImport int32  CALL_CONV_IMP swe_pheno_ut(double tjd_ut, int32 ipl, int32 iflag, double *attr, char *serr);

DllImport double  CALL_CONV_IMP swe_refrac(double inalt, double atpress, double attemp, int32 calc_flag);
DllImport double  CALL_CONV_IMP swe_refrac_extended(double inalt, double geoalt, double atpress, double attemp, double lapse_rate, int32 calc_flag, double *dret);
DllImport void  CALL_CONV_IMP swe_set_lapse_rate(double lapse_rate);

DllImport void  CALL_CONV_IMP swe_azalt(
      double tjd_ut,
      int32 calc_flag,
      double *geopos,
      double atpress,
      double attemp,
      double *xin,
      double *xaz);

DllImport void  CALL_CONV_IMP swe_azalt_rev(
      double tjd_ut,
      int32 calc_flag,
      double *geopos,
      double *xin,
      double *xout);

DllImport int32  CALL_CONV_IMP swe_rise_trans(
               double tjd_ut, int32 ipl, char *starname,
	       int32 epheflag, int32 rsmi,
               double *geopos,
	       double atpress, double attemp,
               double *tret,
               char *serr);

DllImport int32  CALL_CONV_IMP swe_rise_trans_true_hor(
               double tjd_ut, int32 ipl, char *starname,
	       int32 epheflag, int32 rsmi,
               double *geopos,
	       double atpress, double attemp,
	       double horhgt,
               double *tret,
               char *serr);

DllImport int32  CALL_CONV_IMP swe_nod_aps(double tjd_et, int32 ipl, int32 iflag,
                      int32  method,
                      double *xnasc, double *xndsc,
                      double *xperi, double *xaphe,
                      char *serr);

DllImport int32  CALL_CONV_IMP swe_nod_aps_ut(double tjd_ut, int32 ipl, int32 iflag,
                      int32  method,
                      double *xnasc, double *xndsc,
                      double *xperi, double *xaphe,
                      char *serr);

DllImport int32 CALL_CONV_IMP swe_get_orbital_elements(double tjd_et, int32 ipl, int32 iflag, double *dret, char *serr);

DllImport int32 CALL_CONV_IMP swe_orbit_max_min_true_distance(double tjd_et, int32 ipl, int32 iflag, double *dmax, double *dmin, double *dtrue, char *serr);

DllImport centisec  CALL_CONV_IMP swe_csnorm(centisec p);

DllImport centisec  CALL_CONV_IMP swe_difcsn (centisec p1, centisec p2);

DllImport double  CALL_CONV_IMP swe_difdegn (double p1, double p2);

DllImport centisec  CALL_CONV_IMP swe_difcs2n(centisec p1, centisec p2);

DllImport double  CALL_CONV_IMP swe_difdeg2n(double p1, double p2);

DllImport double  CALL_CONV_IMP swe_difdeg2n(double p1, double p2);
DllImport double  CALL_CONV_IMP swe_difrad2n(double p1, double p2);
DllImport double  CALL_CONV_IMP swe_rad_midp(double x1, double x0);
DllImport double  CALL_CONV_IMP swe_deg_midp(double x1, double x0);

DllImport centisec  CALL_CONV_IMP swe_csroundsec(centisec x);

DllImport int32  CALL_CONV_IMP swe_d2l(double x);

DllImport void  CALL_CONV_IMP swe_split_deg(double ddeg, int32 roundflag, int32 *ideg, int32 *imin, int32 *isec, double *dsecfr, int32 *isgn);

DllImport int  CALL_CONV_IMP swe_day_of_week(double jd);

DllImport char * CALL_CONV_IMP swe_cs2timestr(CSEC t, int sep, AS_BOOL suppressZero, char *a);

DllImport char * CALL_CONV_IMP swe_cs2lonlatstr(CSEC t, char pchar, char mchar, char *s);

DllImport char * CALL_CONV_IMP swe_cs2degstr(CSEC t, char *a);

DllImport void CALL_CONV_IMP swe_set_interpolate_nut(AS_BOOL do_interpolate);

#endif
#ifdef __cplusplus
}
#endif
