#define PREC_IAU_1976_CTIES          2.0
#define PREC_IAU_2000_CTIES          2.0

#define PREC_IAU_2006_CTIES          75.0

#define DPSI_DEPS_IAU1980_FILE_EOPC04   "eop_1962_today.txt"
#define DPSI_DEPS_IAU1980_FILE_FINALS   "eop_finals.txt"
#define DPSI_DEPS_IAU1980_TJD0_HORIZONS  2437684.5
#define HORIZONS_TJD0_DPSI_DEPS_IAU1980  2437684.5
#define DPSI_IAU1980_TJD0	(64.284 / 1000.0)
#define DEPS_IAU1980_TJD0	(6.151 / 1000.0)

extern void swi_coortrf(double *xpo, double *xpn, double eps);

extern void swi_coortrf2(double *xpo, double *xpn, double sineps, double coseps);

extern void swi_cartpol(double *x, double *l);

extern void swi_cartpol_sp(double *x, double *l);
extern void swi_polcart_sp(double *l, double *x);

extern void swi_polcart(double *l, double *x);

extern void swi_bias(double *x, double tjd, int32 iflag, AS_BOOL backward);
extern void swi_get_eop_time_range(void);

extern void swi_icrs2fk5(double *x, int32 iflag, AS_BOOL backward);

extern int swi_precess(double *R, double J, int32 iflag, int direction );
extern void swi_precess_speed(double *xx, double t, int32 iflag, int direction);

extern int32 swi_guess_ephe_flag(void);

extern void swi_deflect_light(double *xx, double dt, int32 iflag);
extern void swi_aberr_light(double *xx, double *xe, int32 iflag);
extern int swi_plan_for_osc_elem(int32 iflag, double tjd, double *xx);
extern int swi_trop_ra2sid_lon(double *xin, double *xout, double *xoutr, int32 iflag);
extern int swi_trop_ra2sid_lon_sosy(double *xin, double *xout, int32 iflag);
extern int swi_get_observer(double tjd, int32 iflag,
	AS_BOOL do_save, double *xobs, char *serr);
extern void swi_force_app_pos_etc(void);

extern void swi_check_ecliptic(double tjd, int32 iflag);
extern double swi_epsiln(double J, int32 iflag);
extern void swi_ldp_peps(double J, double *dpre, double *deps);

extern void swi_check_nutation(double tjd, int32 iflag);
extern int swi_nutation(double J, int32 iflag, double *nutlo);
extern void swi_nutate(double *xx, int32 iflag, AS_BOOL backward);

extern void swi_mean_lunar_elements(double tjd,
							 double *node, double *dnode,
							 double *peri, double *dperi);

extern double swi_mod2PI(double x);

extern double swi_echeb(double x, double *coef, int ncf);
extern double swi_edcheb(double x, double *coef, int ncf);

extern void swi_cross_prod(double *a, double *b, double *x);

extern double swi_dot_prod_unit(double *x, double *y);

extern double swi_angnorm(double x);

extern void swi_gen_filename(double tjd, int ipli, char *fname);

extern uint32 swi_crc32(unsigned char *buf, int len);

extern int swi_cutstr(char *s, char *cutlist, char *cpos[], int nmax);
extern char *swi_right_trim(char *s);

extern double swi_kepler(double E, double M, double ecce);

extern char *swi_get_fict_name(int32 ipl, char *s);

extern void swi_FK4_FK5(double *xp, double tjd);

extern char *swi_strcpy(char *to, char *from);
extern char *swi_strncpy(char *to, char *from, size_t n);

extern double swi_deltat_ephe(double tjd_ut, int32 epheflag);

#ifdef TRACE
#  define TRACE_COUNT_MAX         10000
  extern TLS FILE *swi_fp_trace_c;
  extern TLS FILE *swi_fp_trace_out;
  extern TLS int32 swi_trace_count;
  extern void swi_open_trace(char *serr);
  static const char *fname_trace_c = "swetrace.c";
  static const char *fname_trace_out = "swetrace.txt";
#ifdef FORCE_IFLAG
  static const char *fname_force_flg = "force.flg";
#endif
#endif
