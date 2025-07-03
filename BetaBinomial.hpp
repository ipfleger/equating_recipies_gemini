/* BetaBinomial.hpp 
*/ 

#ifndef BETABINOMIAL_HPP
#define BETABINOMIAL_HPP

#include <cstdio> // For FILE
#include <string> // For std::string

// Note: External dependency 'ERutilities.h' is not provided. 
// It is assumed to contain definitions for structs like USTATS, BB_SMOOTH, PDATA, ERAW_RESULTS.
#include "ERutilities.h"

// C-style macros replaced with C++ constexpr constants.
constexpr int RMAXIT = 20;    // maximum number of iterations for computing upper limit
constexpr double KACC = 0.001; // accuracy for computing lower limit in CalcBetaParaLS

// Forward declarations of structs assumed to be in ERutilities.h
struct USTATS;
struct BB_SMOOTH;
struct PDATA;
struct ERAW_RESULTS;

void Wrapper_Smooth_BB(USTATS* x, int nparm, double rel,
                       BB_SMOOTH* s);
void Smooth_BB(int n, int nitems, double* fd, double* mts, 
               int nparm, double rel, BB_SMOOTH* s);
void Print_BB(FILE* fp, const std::string& tt, USTATS* x, BB_SMOOTH* s);
void Wrapper_RB(char design, char method, char smoothing,  
                USTATS* x, USTATS* y, 
                BB_SMOOTH* bbx, BB_SMOOTH* bby, int rep,
                PDATA* inall, ERAW_RESULTS* r);
void Print_RB(FILE* fp, const std::string& tt, PDATA* inall, ERAW_RESULTS* r);

short Beta4Smooth(double* rfreq, BB_SMOOTH* betafit);
short Beta2Smooth(double* rfreq, BB_SMOOTH* betafit);

short CalcBetaPara(int nitems, double* moment, double* nctmoment, double* para);
short EstNegHypGeo(int nitems, double* moment, double* para);
short CalcBetaParaMM(int nitems, double* moment, double* para);
short CalcBetaParaLS(int nitems, double* moment, double* nctmoment, double* para);
short FindUpper(double* para, double* tmoment);
short FindLower(double* para, double* tmoment);
double CalcKurt(double* para);
short Kurtfuncd(double x, int nitems, double* para, double* tmoment, double* nctmoment,
		        double* kurt);
short KurtfuncUpper(double x, int nitems, double* para, double* tmoment, double* nctmoment,
		            double* kurt);
		
void BetaMoments(int n, double k, const double* rmoment, double* tmoment, double* nctmoment);

short ObsDensity(int n, int nsamp, double* beta, double* scounts);
void ObsDenK(int n, double k, double* scounts);
double dotprod(int n, double* v1, double* v2);
void CalcM3(int n, double* beta, double** m3);
short CalcM24(int n, double para, double* m);
void CalcM15(int n, double para, double* m);

double LRChiSqr(int n, const double* rawc, const double* fitc, int* ncat);
double PChiSqr(int n, double minexpcount, const double* rawc,
	             const double* fitc, int* ncat);
	
double CalcLordk(double kr20, int nitems, double* rmoment);

#endif // BETABINOMIAL_HPP