/*
ERutilities.hpp

  Includes:  (a) all structure declarations, and
             (b) prototypes for functions in ERutilities.cpp
*/

#ifndef ERUTILITIES_HPP
#define ERUTILITIES_HPP

// Microsoft-specific pragma to disable certain warnings.
#if defined(_MSC_VER)
#pragma warning(disable:4996)
#pragma warning(disable:4706)
#pragma warning(disable:4715)
#pragma warning(disable:4701)
#pragma warning(disable:4127)
#pragma warning(disable:4244)
#pragma warning(disable:4100)
#pragma warning(disable:4305)
#endif

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>

// Forward declarations to break circular include dependencies.
// The full definitions are in their respective .hpp files.
struct ItemSpec;
struct IRTstControl;
struct RawFitDist;
struct RawTruObsEquiv;
struct IRT_INPUT;
struct BB_SMOOTH;
struct ULL_SMOOTH;
struct BLL_SMOOTH;
struct CS_SMOOTH;

/************************ structure definitions **************************/

struct USTATS{
  char fname[100];
  char id;
  int n;
  double mind;
  double maxd;
  double min;
  double max;
  double inc;
  int ns;
  int *fd;
  double *dbl_fd;
  int *cfd;
  double *rfd;
  double *crfd;
  double *prd;
  double mts[4];
};

struct BSTATS{
  char fname[100];
  int n;
  char id1;
  double mind1;
  double maxd1;
  double min1;
  double max1;
  double inc1;
  int ns1;
  int *fd1;
  double *dbl_fd1;
  int *cfd1;
  double *rfd1;
  double *crfd1;
  double *prd1;
  double mts1[4];
  char id2;
  double mind2;
  double maxd2;
  double min2;
  double max2;
  double inc2;
  int ns2;
  int *fd2;
  double *dbl_fd2;
  int *cfd2;
  double *rfd2;
  double *crfd2;
  double *prd2;
  double mts2[4];
  int **bfd;
  double **dbl_bfd;
  double cov;
  double corr;
  double **bp12;
};

struct PDATA{
  char xfname[100];
  char yfname[100];
  char xyfname[100];
  USTATS *x;
  USTATS *y;
  BSTATS *xv;
  BSTATS *yv;
  BSTATS *xy;
  char design;
  char method;
  char smoothing;
  double w1;
  int anchor;
  double rv1;
  double rv2;
  int nm;
  char **names;
  double min;
  double max;
  double inc;
  int *fdx;
  int n;
  double minp;
  double maxp;
  double incp;
  char nameyct[100];
  double **yct;
  int round;
  int lprss;
  int hprss;
  int nrep;
  int rep;
  BB_SMOOTH *bbx;
  BB_SMOOTH *bby;
  ULL_SMOOTH *ullx;
  ULL_SMOOTH *ully;
  BLL_SMOOTH *bllxv;
  BLL_SMOOTH *bllyv;
  BLL_SMOOTH *bllxy;
  CS_SMOOTH *cs;
  IRT_INPUT *IRT_Input;
};

struct ERAW_RESULTS{
  double msx[4];
  double msy[4];
  double ssx[4];
  double ssy[4];
  double gamma1[4];
  double gamma2[4];
  double a[4];
  double b[4];
  double **eraw;
  double **mts;
  double **fxs;
  double **gys;
};

struct ESS_RESULTS{
  double **essu;
  double **essr;
  double **mtsu;
  double **mtsr;
};

struct BOOT_ERAW_RESULTS{
  double **mn;
  double **sd;
  double *bse;
};

struct BOOT_ESS_RESULTS{
  double **mnu;
  double **sdu;
  double *bseu;
  double **mnr;
  double **sdr;
  double *bser;
};

struct BB_SMOOTH{
  int num_items;
  int num_persons;
  int nparm;
  double rel;
  double lordk;
  double beta[4];
  double rmoments[4];
  double fmoments[4];
  double tmoments[4];
  double lrchisq;
  double pchisq;
  short  momentsfit;
  double *density;
  double *crfd;
  double *prd;
};

struct ULL_SMOOTH{
  int num_persons;
  int ns;
  double min;
  double inc;
  int c;
  double **B_raw;
  double **B;
  double *nct;
  double *mct;
  double *Beta;
  double ap;
  double *n_mts;
  double *m_mts;
  double *n_mts_raw;
  double *m_mts_raw;
  int nit;
  int max_nit;
  int ctype;
  int Btype;
  int scale;
  double crit;
  double lrchisq;
  int nzero;
  double *density;
  double *crfd;
  double *prd;
};

struct BLL_SMOOTH{
  int anchor;
  int num_persons;
  int nsu;
  double minu;
  double incu;
  int nsv;
  double minv;
  double incv;
  int ns;
  int cu;
  int cv;
  int cuv;
  int **cpm;
  int nc;
  double **B_raw;
  double **B;
  double *nct;
  double *mct;
  double *Beta;
  double ap;
  double *n_mts;
  double *m_mts;
  double *n_mts_raw;
  double *m_mts_raw;
  int nit;
  int max_nit;
  int ctype;
  int Btype;
  int scale;
  double crit;
  double lrchisq;
  int nzero;
  int nsx;
  double minx;
  double incx;
  double **bfd;
  double *fd_x;
  double *density_x;
  double *crfd_x;
  double *prd_x;
  double *fd_v;
  double *density_v;
  double *crfd_v;
  double *prd_v;
  double **brfd;
  double *crfd_vector_bfd;
};

struct CS_SMOOTH{
  int ns;
  double s;
  double prlow;
  double prhigh;
  int low;
  int high;
  int nsb;
  double *eeq;
  double *se;
  double *cmat;
  double *eeqs;
  double *inv;
};

/************************ function prototypes ****************************/

int loc(double x, double min, double inc);
int nscores(double max, double min, double inc);
double score(int loc, double min, double inc);
int skipcols(FILE *fp, int k);
int atodouble(FILE *fp, char *str, int begin, int end, double *x);
int atointeger(FILE *fp, char *str, int begin, int end, int *x);
void flushline(FILE *fp);
void runerror(const char* error_text);
void convertFtoW(const char* fname, int nv, int fields[][3], const char* tname);

void cum_rel_freqs(double min, double max, double inc,
                   double *rfd, double *crfd);
double perc_rank(double min, double max, double inc, double *crfd, double x);
double interpolate(double x, int ns, double *f);

int ReadRawGet_moments(const char* fname, int scol,
                       double *moments, double *mind, double *maxd);
int ReadRawGet_mind_maxd(const char* fname, int scol, double *mind, double *maxd);
void ReadFdGet_USTATS(const char* fname, int scol, int fcol, double min,
					  double max, double inc, char id, USTATS *s);
void ReadRawGet_USTATS(const char* fname, int scol, double min,
                       double max, double inc, char id, USTATS *s);
void ReadRawGet_BSTATS(const char* fname, int rows, int cols,
                double rmin, double rmax, double rinc,
                double cmin, double cmax, double cinc,
                char rid, char cid, BSTATS *s);
void ReadSSConvTableForY(const char* nameyct, double minp, double maxp, double inc,
                         double **yct);
int MomentsFromFD(double min, double max, double inc, double *scores,
                  int *fd, double *moments);
void MomentsFromRFD(double min, double max, double inc, double *scores,
                    double *fd, double *moments);

void EquiEquate(int nsy, double miny, double incy,
               double *crfdy, int nsx, double *prdx, double *eraw);
double perc_point(int ns, double min, double inc,
                  double *crfd, double pr);

void Wrapper_ESS(PDATA *inall, ERAW_RESULTS *r, double minp,
                double maxp, double incp, const char *nameyct, int round,
                int lprss, int hprss, ESS_RESULTS *s);
void Equated_ss(double min, double max, double inc, double minp, double maxp,
                double incp, double *eraw, double **yct, int round, int lprss,
                int hprss, double *essu, double *essr);

void Print_ESS(FILE *fp, const std::string& tt, PDATA *inall, ESS_RESULTS *s);
void Print_USTATS(FILE *fp, const std::string& tt, USTATS *s);
void Print_BSTATS(FILE *fp, const std::string& tt, BSTATS *s, int pbfd);
void Print_vector(FILE *fp, const std::string& label, double *vector, int nrows,
                  const std::string& rowhead, const std::string& colhead);
void Print_matrix(FILE *fp, const std::string& label, double **matrix, int nrows,
                  int ncols, const std::string& rowhead, const std::string& colheads);
void Print_file(const char* FileName, FILE *fp);

enum class datatype { CHARACTER, INTEGER, FLOAT, DOUBLE };
void	er_error2(const char* err_source, const char* err_message);
void	release_matrix(datatype mode, void **matrix, long rid_low, long cid_low);
void **	allocate_matrix(datatype mode,long rid_low, long rid_high, long cid_low, long cid_high);

void    er_scale(double vect1[], double scale,double vect2[],int n,int offset);
double  er_dot(double vect1[], double vect2[], int n, int offset);
void    er_daxpy(double vectY[], double scale, double vectX[], int n, int offset);
void    er_r1update(double ** matx,double scale,double vect[],int n,int offset);
void    er_mvmult(double vect2[], double **matx,double vect1[], int n, int offset);

float er_random(long *seed);

void	er_ludcmp(double **a, int n, int *pivot, double *det);
void    er_lubksb(double **a, int n, int pivot[], double b[]);

void	er_sort(float* vector, int left, int right);

void    er_matrix_inverse(int n, double **a);

void	er_qrdcmp(double **a, int nrow, int ncol, double *coeff, double *diag);

double  er_rtsafe(void (*funcd)(double, double *, double *), double x0, double x1, double error);

int     er_lnsrch(double xold[], int n, double gk[], double sk[],  double maxstep,
			     double (*ftn_ptr)(double []), double xnew[]);
void    er_dfpmin(double xold[], int n, double error, int * numiter, double *fvalue,
				  double(*ftn_ptr)(double []), void (*dftn_ptr)(double [], double []));

#endif // ERUTILITIES_HPP