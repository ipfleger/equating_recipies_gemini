#ifndef CUBICSPLINE_HPP
#define CUBICSPLINE_HPP

#include <cmath>
#include <cstdio>
#include <string>

// Note: External dependency 'ERutilities.h' is not provided.
#include "ERutilities.h"

// C-style macro for PI replaced with a C++ constexpr variable.
// C++20 provides std::numbers::pi as an alternative.
constexpr double PI = 3.14159265358979323846;

// Forward declarations of structs assumed to be in ERutilities.h
struct PDATA;
struct ERAW_RESULTS;
struct CS_SMOOTH;

/* functions related to wrapper and print */
void Wrapper_Smooth_CubSpl(char design,
                           PDATA *xtoy, ERAW_RESULTS *r_xtoy,
						   double *se_xtoy, CS_SMOOTH *cs_xtoy,
                           PDATA *ytox, ERAW_RESULTS *r_ytox,
						   double *se_ytox, CS_SMOOTH *cs_ytox,
                           double prlow, double prhigh, double s, int rep,
                           PDATA *inall, ERAW_RESULTS *r);
void Smooth_CubSpl(char design, PDATA *z, ERAW_RESULTS *r,
				   double *se, CS_SMOOTH *cs_xtoy,
				   double prlow, double prhigh, double s,
			       int rep, int ns1, int ns2, int inverse);
void Print_CubSpl(FILE *fp, const std::string& tt, PDATA *xtoy, PDATA *ytox,
                  PDATA *inall, ERAW_RESULTS *r, int parm);
void Print_Coefs_CubSpl(FILE *fp, PDATA *z, char c);

/* functions related to numerical linear algebra */
void	dpdch(double * mx,int dim);
void	dtdmt(double *tdm2,double *tdm1, int rdim1, int cdim1);
void	dtdsm(double *tdm2,double c, double *tdm1, int rdim1, int cdim1);
void	dtdmm(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1, int cdim2);
void	dtdma(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1);
void	subbk(double* b, double *utm, int dim);
void	subfd(double* b, double *utm, int dim);
void	chsol(double *vx, double* vb, double *mpd, int dim);
void	dtdmv(double *v2, double * mtd, double * v1, int rdim, int cdim);

/* functions related to numerical linear/cubic polynomial */
double	linearPoly(double x0, double y0, double x1, double y1, double xvalue);
double	cubicPoly(double xleft, double xright, double ai, double bi, double ci,
		double di, double xvalue);

/* functions related to post-smoothing method using cubic splines */
void	setupT(bool edist, double * mt, int n, double * h);
void	setupQt(bool edist, double * mqt,int n, double * h);
int		sspline(double *x, double *y,double *dyi, int n, double s, double *cmx);
void	postSmooth(double *x, double *y, double *dyi, int num1, double s,
		int xlow, int xhigh, double ky,double *vectX, int num2, double *vectY, double *cmat);

/* functions related to the inverse of the post-smoothing method using cubic splines */
double	inverseCubicPoly(double left, double right, double ai, double bi,
		double ci, double di, double yvalue);
void	inverseSSpline(double *ynodes, double * cmat, int n2,double *vectX, int n1,double *vectY);
void	inversePostSmooth(double *yvalues, double *xvalues, double *dxi, int num1, double s,
		int ylow, int yhigh, double kx,double *vectX, int num2, double *vectY, double *cmat);

#endif // CUBICSPLINE_HPPvvvvv