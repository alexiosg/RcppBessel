// Toplevel Routines called from R :

void zbesh(double *zr, double *zi, double *fnu, int *kode, int *m, int *n, double *cyr, double *cyi, int *nz, int *ierr);
void zbesi(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
void zbesj(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
void zbesk(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
void zbesy(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, double *cwrkr, double *cwrki, int *ierr);
void zairy(double *zr, double *zi, int *id, int *kode, double *air, double *aii, int *nz, int *ierr);
void zbiry(double *zr, double *zi, int *id, int *kode, double *bir, double *bii, int *ierr);

// Auxiliaries

// TODO: declare all  'static void' / 'static <foo>' (not just void or 'int')
void zmlt(double ar, double ai, double br, double bi, double *cr, double *ci);
void zdiv(double ar, double ai, double br, double bi, double *cr, double *ci);
int zsqrt_sub__(double *ar, double *ai, double *br, double *bi);
int zexp_sub__(double *ar, double *ai, double *br, double *bi);
int zlog_sub__(double *ar, double *ai, double *br, double *bi, int *ierr);
double zabs(double zr, double zi);
int zbknu(double *zr, double *zi, double *fnu, int kode, int n, int verbose,
	  double *yr, double *yi, double tol, double elim, double alim);
int zkscl_(double *zrr, double *zri, double *fnu, int *n, double *yr, double *yi, int *nz, double *rzr, double *rzi, double *ascle, double *tol, double *elim);
int zshch_(double *zr, double *zi, double *cshr, double *cshi, double *cchr, double *cchi);
void zrati(double zr, double zi, double fnu, int n, double tol,
	   double *cyr, double *cyi);
void zs1s2_(double *zrr, double *zri, double *s1r, double *s1i, double *s2r, double *s2i, int *nz, double *ascle, double *alim, int *iuf);
int zbunk_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zmlri(double *zr, double *zi, double *fnu, int kode, int n, double *yr, double *yi, double tol);
void zwrsk_(double *zrr, double *zri, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *cwr, double *cwi, double *tol, double *elim, double *alim);
void zseri_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zasyi(double *zr, double *zi, double *fnu, int kode, int n, double *yr, double *yi,
	  double rl, double tol, double elim, double alim);
int zuoik(double *zr, double *zi, double *fnu, int kode, int ikflg, int n, double *yr, double *yi, double tol, double elim, double alim);
void zacon_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz,
	    int verbose, double rl, double fnul, double tol, double elim, double alim);
int zbinu(double *zr, double *zi, double *fnu, int kode, int n, double *cyr, double *cyi,
	  double rl, double fnul, double tol, double elim, double alim);
double dgamln(double z, int *ierr);

void zacai_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *rl, double *tol, double *elim, double *alim);
int zuchk(double yr, double yi, double ascle, double tol);
int zunik_(double *zrr, double *zri, double *fnu, int *ikflg, int *ipmtr, double *tol, int *init, double *phir, double *phii, double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i, double *sumr, double *sumi, double *cwrkr, double *cwrki);
int zunhj_(double *zr, double *zi, double *fnu, int *ipmtr, double *tol, double *phir, double *phii, double *argr, double *argi, double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i, double *asumr, double *asumi, double *bsumr, double *bsumi);
void zunk1_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
void zunk2_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zbuni(double *zr, double *zi, double *fnu,
	  int kode, int n, int nui,
	  double *yr, double *yi, int *nlast,
	  double fnul, double tol, double elim, double alim);
void zuni1_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, int *nlast, double *fnul, double *tol, double *elim, double *alim);
void zuni2_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, int *nlast, double *fnul, double *tol, double *elim, double *alim);

/*:ref: d1mach_ 7 1 4 */
/*:ref: i1mach_ 4 1 4 */
