#define SE_VERSION      "2.10.03"

#define J2000           2451545.0
#define B1950           2433282.42345905
#define J1900           2415020.0
#define B1850           2396758.2035810

#define MPC_CERES       1
#define MPC_PALLAS      2
#define MPC_JUNO        3
#define MPC_VESTA       4
#define MPC_CHIRON      2060
#define MPC_PHOLUS      5145

#define SE_NAME_SUN             "Sun"
#define SE_NAME_MOON            "Moon"
#define SE_NAME_MERCURY         "Mercury"
#define SE_NAME_VENUS           "Venus"
#define SE_NAME_MARS            "Mars"
#define SE_NAME_JUPITER         "Jupiter"
#define SE_NAME_SATURN          "Saturn"
#define SE_NAME_URANUS          "Uranus"
#define SE_NAME_NEPTUNE         "Neptune"
#define SE_NAME_PLUTO           "Pluto"
#define SE_NAME_MEAN_NODE       "mean Node"
#define SE_NAME_TRUE_NODE       "true Node"
#define SE_NAME_MEAN_APOG       "mean Apogee"
#define SE_NAME_OSCU_APOG       "osc. Apogee"
#define SE_NAME_INTP_APOG       "intp. Apogee"
#define SE_NAME_INTP_PERG       "intp. Perigee"
#define SE_NAME_EARTH           "Earth"
#define SE_NAME_CERES           "Ceres"
#define SE_NAME_PALLAS          "Pallas"
#define SE_NAME_JUNO            "Juno"
#define SE_NAME_VESTA           "Vesta"
#define SE_NAME_CHIRON          "Chiron"
#define SE_NAME_PHOLUS          "Pholus"

#define SE_NAME_CUPIDO          "Cupido"
#define SE_NAME_HADES           "Hades"
#define SE_NAME_ZEUS            "Zeus"
#define SE_NAME_KRONOS          "Kronos"
#define SE_NAME_APOLLON         "Apollon"
#define SE_NAME_ADMETOS         "Admetos"
#define SE_NAME_VULKANUS        "Vulkanus"
#define SE_NAME_POSEIDON        "Poseidon"
#define SE_NAME_ISIS            "Isis"
#define SE_NAME_NIBIRU          "Nibiru"
#define SE_NAME_HARRINGTON      "Harrington"
#define SE_NAME_NEPTUNE_LEVERRIER       "Leverrier"
#define SE_NAME_NEPTUNE_ADAMS   "Adams"
#define SE_NAME_PLUTO_LOWELL    "Lowell"
#define SE_NAME_PLUTO_PICKERING "Pickering"
#define SE_NAME_VULCAN          "Vulcan"
#define SE_NAME_WHITE_MOON      "White Moon"

#define PI              M_PI
#define TWOPI           (2.0 * PI)

#define ENDMARK         -99

#define SEI_EPSILON     -2
#define SEI_NUTATION    -1
#define SEI_EMB		0
#define SEI_EARTH	0
#define SEI_SUN  	0
#define	SEI_MOON	1
#define	SEI_MERCURY	2
#define	SEI_VENUS	3
#define	SEI_MARS	4
#define	SEI_JUPITER	5
#define	SEI_SATURN	6
#define	SEI_URANUS	7
#define	SEI_NEPTUNE	8
#define	SEI_PLUTO	9
#define	SEI_SUNBARY	10
#define	SEI_ANYBODY	11
#define	SEI_CHIRON	12
#define	SEI_PHOLUS	13
#define	SEI_CERES	14
#define	SEI_PALLAS	15
#define	SEI_JUNO	16
#define	SEI_VESTA	17

#define SEI_NPLANETS    18

#define SEI_MEAN_NODE   0
#define SEI_TRUE_NODE   1
#define SEI_MEAN_APOG   2
#define SEI_OSCU_APOG   3
#define SEI_INTP_APOG   4
#define SEI_INTP_PERG   5

#define SEI_NNODE_ETC    6

#define SEI_FLG_HELIO   1
#define SEI_FLG_ROTATE  2
#define SEI_FLG_ELLIPSE 4
#define SEI_FLG_EMBHEL  8

#define SEI_FILE_PLANET	  0
#define SEI_FILE_MOON	  1
#define SEI_FILE_MAIN_AST 2
#define SEI_FILE_ANY_AST  3
#define SEI_FILE_FIXSTAR  4
#define SEI_FILE_PLMOON   5

#if 0
#define SEI_FILE_TEST_ENDIAN     (97L * 65536L + 98L * 256L + 99L)
#endif
#define SEI_FILE_TEST_ENDIAN     (0x616263L)
#define SEI_FILE_BIGENDIAN	0
#define SEI_FILE_NOREORD	0
#define SEI_FILE_LITENDIAN	1
#define SEI_FILE_REORD  	2

#define SEI_FILE_NMAXPLAN	50
#define SEI_FILE_EFPOSBEGIN      500

#define SE_FILE_SUFFIX	"se1"

#define SEI_NEPHFILES   7
#define SEI_CURR_FPOS   -1
#define SEI_NMODELS 8

#define SEI_ECL_GEOALT_MAX   25000.0
#define SEI_ECL_GEOALT_MIN   (-500.0)

#define CHIRON_START    1967601.5
#define CHIRON_END      3419437.5

#define PHOLUS_START    640648.5
#define PHOLUS_END      4390617.5

#define MOSHPLEPH_START	 625000.5
#define MOSHPLEPH_END  	2818000.5
#define MOSHLUEPH_START	 625000.5
#define MOSHLUEPH_END  	2818000.5

#define MOSHNDEPH_START	-3100015.5
#define MOSHNDEPH_END  	8000016.5

#define JPL_DE431_START -3027215.5
#define JPL_DE431_END    7930192.5

#if FALSE
#define JPLEPH_START	 625307.5
#define JPLEPH_END	2816848.5
#define SWIEPH_START	 625614.927151
#define SWIEPH_END	2813641.5
#define ALLEPH_START	MOSHPLEPH_START
#define ALLEPH_END	MOSHPLEPH_END
#define BEG_YEAR       (-3000)
#define END_YEAR       3000
#endif

#define MAXORD          40

#define NCTIES         6.0

#define OK (0)
#define ERR (-1)
#define NOT_AVAILABLE (-2)
#define BEYOND_EPH_LIMITS (-3)

#define J_TO_J2000   	1
#define J2000_TO_J   	-1

#define MOON_MEAN_DIST  384400000.0
#define MOON_MEAN_INCL  5.1453964
#define MOON_MEAN_ECC   0.054900489

#define SUN_EARTH_MRAT  332946.050895
#define EARTH_MOON_MRAT (1 / 0.0123000383)
#if 0
#define EARTH_MOON_MRAT 81.30056907419062
#endif
#if 0
#define EARTH_MOON_MRAT 81.30056
#endif

#define AUNIT       	1.49597870700e+11
#define CLIGHT       	2.99792458e+8
#if 0
#define HELGRAVCONST    1.32712438e+20
#endif
#define HELGRAVCONST    1.32712440017987e+20
#define GEOGCONST       3.98600448e+14
#define KGAUSS		0.01720209895
#define SUN_RADIUS      (959.63 / 3600 * DEGTORAD)
#define EARTH_RADIUS	6378136.6

#define EARTH_OBLATENESS (1.0/ 298.25642)
#define EARTH_ROT_SPEED (7.2921151467e-5 * 86400)

#define LIGHTTIME_AUNIT  (499.0047838362/3600.0/24.0)
#define PARSEC_TO_AUNIT  206264.8062471

#define SSY_PLANE_NODE_E2000    (107.582569 * DEGTORAD)

#define SSY_PLANE_NODE          (107.58883388 * DEGTORAD)

#define SSY_PLANE_INCL          (1.578701 * DEGTORAD)

#define KM_S_TO_AU_CTY	 21.095
#define MOON_SPEED_INTV  0.00005
#define PLAN_SPEED_INTV  0.0001
#define MEAN_NODE_SPEED_INTV  0.001
#define NODE_CALC_INTV  0.0001
#define NODE_CALC_INTV_MOSH   0.1
#define NUT_SPEED_INTV   0.0001
#define DEFL_SPEED_INTV  0.0000005

#define SE_LAPSE_RATE        0.0065

#define square_sum(x)   (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define dot_prod(x,y)   (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])

#define PNOINT2JPL {J_EARTH, J_MOON, J_MERCURY, J_VENUS, J_MARS, J_JUPITER, J_SATURN, J_URANUS, J_NEPTUNE, J_PLUTO, J_SUN, }

#define NDIAM  (SE_VESTA + 1)
static const double pla_diam[NDIAM] = {1392000000.0,
                           3475000.0,
                           2439400.0 * 2,
                           6051800.0 * 2,
                           3389500.0 * 2,
                          69911000.0 * 2,
                          58232000.0 * 2,
                          25362000.0 * 2,
                          24622000.0 * 2,
                           1188300.0 * 2,
                           0, 0, 0, 0,
                           6371008.4 * 2,
                            271370.0,
                            290000.0,
                            939400.0,
                            545000.0,
                            246596.0,
                            525400.0,
                        };

struct aya_init {double t0;
                 double ayan_t0;
		 AS_BOOL t0_is_UT;
		 int prec_offset;};
static const struct aya_init ayanamsa[SE_NSIDM_PREDEF] = {

{2433282.42346, 24.042044444, FALSE, SEMOD_PREC_NEWCOMB},

{2435553.5, 23.250182778 - 0.004658035, FALSE, SEMOD_PREC_IAU_1976},

{1721057.5, 0, TRUE, 0},

{J1900, 360 - 338.98556, FALSE, SEMOD_PREC_NEWCOMB},

{J1900, 360 - 341.33904, FALSE, -1},

{J1900, 360 - 337.636111, FALSE, SEMOD_PREC_NEWCOMB},

{J1900, 360 - 333.0369024, FALSE, 0},

{J1900, 360 - 338.917778, FALSE, -1},

{J1900, 360 - 338.634444, FALSE, -1},

{1684532.5, -5.66667, TRUE, -1},
{1684532.5, -4.26667, TRUE, -1},
{1684532.5, -3.41667, TRUE, -1},

{1684532.5, -4.46667, TRUE, -1},

{1673941, -5.079167, TRUE, -1},

{1684532.5, -4.44138598, TRUE, 0},

{1674484.0, -9.33333, TRUE, -1},

{1927135.8747793, 0, TRUE, -1},

{0, 0, FALSE, 0},

{J2000, 0, FALSE, 0},

{J1900, 0, FALSE, 0},

{B1950, 0, FALSE, 0},

{1903396.8128654, 0, TRUE, 0},

{1903396.8128654,-0.21463395, TRUE, 0},

{1903396.7895321, 0, TRUE, 0},

{1903396.7895321,-0.23763238, TRUE, 0},

{1903396.8128654,-0.79167046, TRUE, 0},

{1903396.8128654, 2.11070444, TRUE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{2451079.734892000, 30, FALSE, 0},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{1911797.740782065, 0, TRUE, 0},

{1721057.5, -3.2, TRUE, -1},

{0, 0, FALSE, 0},

{0, 0, FALSE, 0},

{2451544.5, 25.0, TRUE, 0},

{1775845.5, -2.9422, TRUE, -1},

{J1900, 22.44597222, FALSE, SEMOD_PREC_NEWCOMB},

{1825235.2458513028, 0.0, FALSE, 0},

{1827424.752255678, 0.0, FALSE, 0},

{2435553.5, 23.25 - 0.00464207, FALSE, SEMOD_PREC_NEWCOMB},

    };

#define PLAN_DATA struct plan_data

struct epsilon {
  double teps, eps, seps, ceps;
};

struct plan_data {

  int ibdy;
  int32 iflg;

  int ncoe;

  int32 lndx0;
  int32 nndx;
  double tfstart;
  double tfend;
  double dseg;

  double telem;
  double prot;
  double qrot;
  double dprot;
  double dqrot;
  double rmax;

  double peri;
  double dperi;
  double *refep;

  double tseg0, tseg1;
  double *segp;

  int neval;

  double teval;
  int32 iephe;
  double x[6];
  int32 xflgs;
  double xreturn[24];

};

#define STR             4.8481368110953599359e-6

extern int swi_mean_node(double jd, double *x, char *serr);
extern int swi_mean_apog(double jd, double *x, char *serr);
extern int swi_moshmoon(double tjd, AS_BOOL do_save, double *xpm, char *serr) ;
extern int swi_moshmoon2(double jd, double *x);
extern int swi_intp_apsides(double J, double *pol, int ipli);

extern int swi_moshplan(double tjd, int ipli, AS_BOOL do_save, double *xpret, double *xeret, char *serr);
extern int swi_moshplan2(double J, int iplm, double *pobj);
extern int swi_osc_el_plan(double tjd, double *xp, int ipl, int ipli, double *xearth, double *xsun, char *serr);
extern FILE *swi_fopen(int ifno, char *fname, char *ephepath, char *serr);
extern int32 swi_init_swed_if_start(void);
extern int32 swi_set_tid_acc(double tjd_ut, int32 iflag, int32 denum, char *serr);
extern int32 swi_get_tid_acc(double tjd_ut, int32 iflag, int32 denum, int32 *denumret, double *tid_acc, char *serr);

extern int32 swi_get_ayanamsa_ex(double tjd_et, int32 iflag, double *daya, char *serr);
extern int32 swi_get_ayanamsa_ex_ut(double tjd_ut, int32 iflag, double *daya, char *serr);
extern int32 swi_get_ayanamsa_with_speed(double tjd_et, int32 iflag, double *daya, char *serr);

extern double swi_armc_to_mc(double armc, double eps);

extern int32 swi_get_denum(int32 ipli, int32 iflag);

struct nut {
  double tnut;
  double nutlo[2];
  double snut, cnut;
  double matrix[3][3];
};

struct plantbl {
  char max_harmonic[9];
  char max_power_of_t;
  signed char *arg_tbl;
  double *lon_tbl;
  double *lat_tbl;
  double *rad_tbl;
  double distance;
};

struct file_data {
  char fnam[AS_MAXCH];
  int fversion;
  char astnam[50];
  int32 sweph_denum;

  FILE *fptr;
  double tfstart;
  double tfend;
  int32 iflg;
  short npl;
  int ipl[SEI_FILE_NMAXPLAN];
};

struct gen_const {
 double clight,
	aunit,
	helgravconst,
	ratme,
	sunradius;
};

struct save_positions {
  int ipl;
  double tsave;
  int32 iflgsave;

  double xsaves[24];
};

struct node_data {

  double teval;
  int32 iephe;
  double x[6];
  int32 xflgs;
  double xreturn[24];

};

struct topo_data {
  double geolon, geolat, geoalt;
  double teval;
  double tjd_ut;
  double xobs[6];
};

struct sid_data {
  int32 sid_mode;
  double ayan_t0;
  double t0;
  AS_BOOL t0_is_UT;
};

#define SWI_STAR_LENGTH 40
struct fixed_star {
  char skey[SWI_STAR_LENGTH + 2];
  char starname[SWI_STAR_LENGTH + 1];
  char starbayer[SWI_STAR_LENGTH + 1];
  char starno[10];
  double epoch, ra, de, ramot, demot, radvel, parall, mag;
};

#define SWE_DATA_DPSI_DEPS  36525

struct interpol {
  double tjd_nut0, tjd_nut2;
  double nut_dpsi0, nut_dpsi1, nut_dpsi2;
  double nut_deps0, nut_deps1, nut_deps2;
};

struct swe_data {
  AS_BOOL ephe_path_is_set;
  AS_BOOL jpl_file_is_open;
  FILE *fixfp;
  char ephepath[AS_MAXCH];
  char jplfnam[AS_MAXCH];
  int32 jpldenum;
  int32 last_epheflag;
  AS_BOOL geopos_is_set;
  AS_BOOL ayana_is_set;
  AS_BOOL is_old_starfile;
  double eop_tjd_beg;
  double eop_tjd_beg_horizons;
  double eop_tjd_end;
  double eop_tjd_end_add;
  int eop_dpsi_loaded;
  double tid_acc;
  AS_BOOL is_tid_acc_manual;
  AS_BOOL init_dt_done;
  AS_BOOL swed_is_initialised;
  AS_BOOL delta_t_userdef_is_set;
  double delta_t_userdef;
  double ast_G;
  double ast_H;
  double ast_diam;
  char astelem[AS_MAXCH * 10];
  int i_saved_planet_name;
  char saved_planet_name[80];

  double *dpsi;
  double *deps;
  int32 timeout;
  int32 astro_models[SEI_NMODELS];
  AS_BOOL do_interpolate_nut;
  struct interpol interpol;
  struct file_data fidat[SEI_NEPHFILES];
  struct gen_const gcdat;
  struct plan_data pldat[SEI_NPLANETS];
#if 0
  struct node_data nddat[SEI_NNODE_ETC];
#else
  struct plan_data nddat[SEI_NNODE_ETC];
#endif
  struct save_positions savedat[SE_NPLANETS+1];
  struct epsilon oec;
  struct epsilon oec2000;
  struct nut nut;
  struct nut nut2000;
  struct nut nutv;
  struct topo_data topd;
  struct sid_data sidd;
  AS_BOOL n_fixstars_real;
  AS_BOOL n_fixstars_named;
  AS_BOOL n_fixstars_records;
  struct fixed_star *fixed_stars;
};

extern TLS struct swe_data swed;
