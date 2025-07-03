/*
  CLL_Equate.hpp
*/

#ifndef CLL_EQUATE_HPP
#define CLL_EQUATE_HPP

#include <cstdio>
#include <string>

// Note: External dependency 'ERutilities.h' is not provided.
#include "ERutilities.h"

// Forward declarations of structs assumed to be in ERutilities.h
struct BLL_SMOOTH;
struct USTATS;
struct ULL_SMOOTH;
struct PDATA;
struct ERAW_RESULTS;
struct BSTATS;

// C-style function pointers. In modern C++, std::function could be an alternative.
typedef double (*PTR_FTN2)(int, double*, double );
typedef double (*PTR_FTN3)(int, double*, double, int );
typedef double (*PTR_FTN4)(double, double, int, double *, double);
typedef double (*PTR_FTN5)(BLL_SMOOTH *, double, double);


double CLLEGPdf(double min, double max, int npara, double *para,
					double x, double nc);

double CLLEGCdf(double min, double max, int npara, double *para,
					double x, double nc);

double GaussianQuadrature16(PTR_FTN2 func, double a, double b, int npara,
						  double *para);

double GaussianQuadrature32(PTR_FTN2 func, double a, double b, int npara,
						  double *para);

double GaussianQuadrature64(PTR_FTN2 func, double a, double b, int npara,
						  double *para);

double GaussianQuadrature64i(PTR_FTN3 func, double a, double b, int npara,
						  double *para,  int i);

double ExpPolynomial(int npara, double *para, double x);

double ExpPolynomialxi(int npara, double *para, double x, int i);

void CalcLLContinuMoments(PTR_FTN4 pdf, double a, double b, int npara,
						  double *para, double *moments);
void CLLEquateEG(double minx, double maxx, int nparax, double *parax,
			   double miny, double maxy, int nparay, double *paray,
			   int nDistCatx, double *scoresx, double *Equatedx);

double CLLInverseCdf(double min, double max, int npara, double *para,
							  double cdf, double nc);

double BivExpPolynomial(BLL_SMOOTH *bivar, double x, double y);

double BivGaussianQuadrature64(PTR_FTN5 func, BLL_SMOOTH *bivar, double ax,
							   double bx, double ay, double by);

double BivGaussianQuadrature32(PTR_FTN5 func, BLL_SMOOTH *bivar, double ax, double bx,
							   double ay, double by);

double CLLMargYPdf(BLL_SMOOTH *bivar, double y, double nc);

double CLLMargXPdf(BLL_SMOOTH *bivar, double x, double nc);

double CLLBivPdf(BLL_SMOOTH *bivar, double x, double y);

double CLLBivCdf(BLL_SMOOTH *bivar, double x, double y, double nc);

double CLLNEATPSMargPdf(BLL_SMOOTH *bivar1, BLL_SMOOTH *bivar2, double *fv,
						double wts, double x);

double CLLNEATPSMargCdf(BLL_SMOOTH *bivar1, BLL_SMOOTH *bivar2, double wts,
								 double x);

double CLLNEATPSInverseCdf(BLL_SMOOTH *bivar1, BLL_SMOOTH *bivar2,
										double wts, double cdf);

int CLLEquateNEATPS(BLL_SMOOTH *bivar1, BLL_SMOOTH *bivar2, double wts,
						   double *Equatedx);

double CLLMargInverseXCdf(BLL_SMOOTH *bivar1, double xcdf, double nc);

double CLLMargInverseYCdf(BLL_SMOOTH *bivar1, double ycdf, double nc);

double CLLMargInverseCBYCdf(BLL_SMOOTH *bivar1, BLL_SMOOTH *bivar2, double wtsy, double ycdf, double nc1, double nc2);

int CLLEquateNEATChn(BLL_SMOOTH *bivar1, BLL_SMOOTH *bivar2, double *Equatedx);

int CLLEquateSG(BLL_SMOOTH *bivar, double *Equatedx);

void CLLSEEEG(long npx, int nparax, double *parax, int nCatx, double *scoresx,
			  long npy, int nparay, double *paray,	int nCaty, double *scoresy, double *SEE);

void Wrapper_RC(char design, char method, char smoothing,
               USTATS *x, USTATS *y,
               ULL_SMOOTH *ullx, ULL_SMOOTH *ully, int rep,
               PDATA *inall, ERAW_RESULTS *r);

void Print_RC(FILE *fp, const std::string& tt, PDATA *inall, ERAW_RESULTS *r);

void Wrapper_SC(char design, char method, char smoothing,  BSTATS *xy,
               BLL_SMOOTH *bllxy, int rep, PDATA *inall,
			   ERAW_RESULTS *r);

void Print_SC(FILE *fp, const std::string& tt, PDATA *inall, ERAW_RESULTS *r);

void Wrapper_CC(char design, char method, char smoothing,
                double w1, int anchor, double rv1, double rv2,
                BSTATS *xv, BSTATS *yv,
				BLL_SMOOTH *bllxv, BLL_SMOOTH *bllyv,
			    int rep, PDATA *inall, ERAW_RESULTS *r);

void Print_CC(FILE *fp, const std::string& tt, PDATA *inall, ERAW_RESULTS *r);

#endif // CLL_EQUATE_HPP