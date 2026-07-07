struct houses {
	  double cusp[37];
	  double cusp_speed[37];
	  double ac;
	  double ac_speed;
	  double mc;
	  double mc_speed;
	  double armc_speed;
	  double vertex;
	  double vertex_speed;
	  double equasc;
	  double equasc_speed;
	  double coasc1;
	  double coasc1_speed;
	  double coasc2;
	  double coasc2_speed;
	  double polasc;
	  double polasc_speed;
	  double sundec;
	  AS_BOOL do_speed;
	  AS_BOOL do_hspeed;
	  AS_BOOL do_interpol;
	  char serr[AS_MAXCH];
	};

#define HOUSES 	struct houses
#define VERY_SMALL	1E-10

#define degtocs(x)    (d2l((x) * DEG))
#define cstodeg(x)    (double)((x) * CS2DEG)

#define sind(x) sin((x) * DEGTORAD)
#define cosd(x) cos((x) * DEGTORAD)
#define tand(x) tan((x) * DEGTORAD)
#define asind(x) (asin(x) * RADTODEG)
#define acosd(x) (acos(x) * RADTODEG)
#define atand(x) (atan(x) * RADTODEG)
#define atan2d(y, x) (atan2(y, x) * RADTODEG)
