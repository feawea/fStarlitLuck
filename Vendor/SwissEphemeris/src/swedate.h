#ifdef __cplusplus
extern "C" {
#endif

#ifndef _SWEDLL_H
extern EXP32 int swe_date_conversion (
	int y , int m , int d ,
     	double utime,
     	char c,
 	double *tgmt);

extern EXP32 double *swe_julday(
	int year, int month, int day, double hour,
	int gregflag);

extern EXP32 void swe_revjul (
	double jd,
	int gregflag,
     	int *jyear, int *jmon, int *jday, double *jut);
#endif
#ifdef __cplusplus
}
#endif
