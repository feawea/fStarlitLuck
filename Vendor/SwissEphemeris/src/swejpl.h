#include "sweodef.h"

#define J_MERCURY	0
#define J_VENUS		1
#define J_EARTH		2
#define J_MARS		3
#define J_JUPITER	4
#define J_SATURN	5
#define J_URANUS	6
#define J_NEPTUNE	7
#define J_PLUTO		8
#define J_MOON		9
#define J_SUN		10
#define J_SBARY		11
#define J_EMB		12
#define J_NUT		13
#define J_LIB		14

extern int swi_pleph(double et, int ntarg, int ncent, double *rrd, char *serr);

extern void swi_close_jpl_file(void);

extern int swi_open_jpl_file(double *ss, char *fname, char *fpath, char *serr);

extern int32 swi_get_jpl_denum(void);

extern void swi_IERS_FK5(double *xin, double *xout, int dir);
