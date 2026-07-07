#define TTINY 0.000001
#define TENTHOFSEC  (1.0 / 86400.0)
#define HUNDTHOFSEC  (1.0 / 864000.0)

#define NEVENTMAX 100000
#define EVENT struct event
struct event {
  double tjd;
  int32  evtype;
  int32  ipla;
  int32  iplb;
  char stnama[40];
  char stnamb[40];
  int32  iasp;
  int32  bpind;
  double  dasp;
  double  dang;
  int32  isign;
  int32  backward;
  double dorb;
  double dret;
};

struct aspdat {
  double tjd;
  int32 iasp;
  double tjd_pre;
  double tjd_post;
  int32 fpos_tjd_post;
};

#define CTYP_MASPECTS	0
#define CTYP_INGRESSES	1
#define CTYP_TRANSITS	2
#define CTYP_VOC	3

#define SPLAN_INGRESS   "0123456789mtAFD"
#define SPLAN_ASPECTS   "0123456789mtAFD,a[136199],f[Gal]"

#define SASP_ASPECTS   "1234567"

#define VOC		struct voc

VOC {
  double tvoc;
  double tingr;
  double tingr0;

  char casp, cpl;
  int isign_voc;
  int isign_ingr;
  int isign_ingr0;

};

#define INGRESS	struct ingress

INGRESS {
  int ipl;
  double tingr;
  int isign;
  int direction;
  int ino;
};
