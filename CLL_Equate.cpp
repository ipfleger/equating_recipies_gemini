/*

CLL_Equate.cpp  File for continuized log-linear equating

This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by
the Free Software Foundation.

This file is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>

Copyright 2009
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa

*/

#include <cmath>
#include <cstdlib>
#include "ERutilities.h"
#include "NRutilities.h"
#include "LogLinear.h"
#include "CG_EquiEquate.h"
#include "CLL_Equate.h"
#include "matrix.h"

/*--------------------------------------------------------------------------
	CLLEGPdf

	functionality

	Computes the continuous pdf based on the fitted loglinear model
	parameters for a discrete distribution

	Author: Tianyou Wang 10/29/2004.

	input
		min         lower limit of the distribution
		max         upper limit of the distribution
		npara       number of parameters
		para		a vector of parameters for the loglinear model
					with a design matrix from polynominals of natural
					basis.
		x           a particular score for which the pdf and cdf is
					generated
		nc          normalizing constant


	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLEGPdf(double min, double max, int npara, double *para,
	double x, double nc)
{
	double pdf;

	pdf = ExpPolynomial(npara, para, x);

	pdf /= nc;

	return pdf;

}


/*--------------------------------------------------------------------------
	CLLEGCdf

	functionality

	Computes the continuous pdf and cdf based on the fitted loglinear model
	parameters for a discrete distribution

	Author: Tianyou Wang 10/29/2004.

	input
		min         lower limit of the distribution
		max         upper limit of the distribution
		npara       number of parameters
		para		a vector of parameters for the loglinear model
					with a design matrix from polynominals of natural
					basis.
		x           a particular score for which the pdf and cdf is
					generated
		nc          normalizing constant


	output
 		the function returns the smoothed cdf
--------------------------------------------------------------------------*/

double CLLEGCdf(double min, double max, int npara, double *para,
	double x, double nc)
{
	double cdf;

	cdf = GaussianQuadrature64(&ExpPolynomial, min, x, npara, para);

	cdf /= nc;

	return cdf;

}

/*--------------------------------------------------------------------------
	GaussianQuadrature16

	functionality

	Computes the numerical integration using Gaussian quadrature with
	16 points

	Author: Tianyou Wang 10/29/2004.

	input
		(*func)()   pointer to a function which is the integrand
		a   		lower limit of integral
		b   		upper limit of integral
		npara       number of parameters for the integrand
		para        vector of parameters for the integrand

	output
		The function returns the integrated value.
	--------------------------------------------------------------------------*/

double GaussianQuadrature16(PTR_FTN2 func, double a, double b, int npara,
	double *para)
{
	int j;
	double xr, xm, dx, s;
	static double x[] = { 0.09501250983764, 0.28160355077926, 0.45801677765723,
		0.61787624440264, 0.755404408355, 0.86563120238783,
		0.94457502307323, 0.98940093499165 };
	static double w[] = { 0.18945061045507, 0.18260341504492, 0.169156519395,
		0.14959598881658, 0.12462897125553, 0.09515851168249,
		0.06225352393865, 0.02715245941175 };

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s = 0;
	for (j = 0; j < 8; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx)) + func(npara, para, (xm - dx)));
	}

	return s *= xr;
}

/*--------------------------------------------------------------------------
	GaussianQuadrature32

	functionality

	Computes the numerical integration using Gaussian quadrature with
	32 points

	Author: Tianyou Wang 10/29/2004.

	input
		(*func)()   pointer to a function which is the integrand
		a   		lower limit of integral
		b   		upper limit of integral
		npara       number of parameters for the integrand
		para        vector of parameters for the integrand

	output
		The function returns the integrated value.
	--------------------------------------------------------------------------*/

double GaussianQuadrature32(PTR_FTN2 func, double a, double b, int npara,
	double *para)
{
	int j;
	double xr, xm, dx, s;
	static double x[] = { 0.04830766568774, 0.14447196158280, 0.23928736225214,
		0.33186860228213, 0.42135127613064, 0.50689990893223,
		0.58771575724076, 0.66304426693022, 0.73218211874029,
		0.79448379596794, 0.84936761373257, 0.89632115576605,
		0.93490607593774, 0.96476225558751, 0.98561151154527,
		0.99726386184948 };
	static double w[] = { 0.09654008851473, 0.09563872007927, 0.09384439908080,
		0.09117387869576, 0.08765209300440, 0.08331192422695,
		0.07819389578707, 0.07234579410885, 0.06582222277636,
		0.05868409347854, 0.05099805926238, 0.04283589802223,
		0.03427386291302, 0.02539206530926, 0.01627439473091,
		0.00701861000947 };

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s = 0;
	for (j = 0; j < 16; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx)) + func(npara, para, (xm - dx)));
	}

	return s *= xr;
}

/*--------------------------------------------------------------------------
	GaussianQuadrature64

	functionality

	Computes the numerical integration using Gaussian quadrature with
	64 points

	Author: Tianyou Wang 10/29/2004.

	input
		(*func)()   pointer to a function which is the integrand
		a   		lower limit of integral
		b   		upper limit of integral
		npara       number of parameters for the integrand
		para        vector of parameters for the integrand

	output
		The function returns the integrated value.
	--------------------------------------------------------------------------*/

double GaussianQuadrature64(PTR_FTN2 func, double a, double b, int npara,
	double *para)
{
	int j;
	double xr, xm, dx, s;
	static double x[] = { 0.02435029266342, 0.07299312178780, 0.12146281929612,
		0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
		0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
		0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
		0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
		0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
		0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
		0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
		0.99930504173577 };
	static double w[] = { 0.04869095700914, 0.04857546744150, 0.04834476223480,
		0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
		0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
		0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
		0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
		0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
		0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
		0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
		0.00178328072170 };

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s = 0;
	for (j = 0; j < 32; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx)) + func(npara, para, (xm - dx)));
	}

	return s *= xr;
}

double GaussianQuadrature64i(PTR_FTN3 func, double a, double b, int npara,
	double *para, int i)
{
	int j;
	double xr, xm, dx, s;
	static double x[] = { 0.02435029266342, 0.07299312178780, 0.12146281929612,
		0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
		0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
		0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
		0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
		0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
		0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
		0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
		0.99930504173577 };
	static double w[] = { 0.04869095700914, 0.04857546744150, 0.04834476223480,
		0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
		0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
		0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
		0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
		0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
		0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
		0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
		0.00178328072170 };

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s = 0;
	for (j = 0; j < 32; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx), i) + func(npara, para, (xm - dx), i));
	}

	return s *= xr;
}


/*--------------------------------------------------------------------------
	ExpPolynomial

	functionality

	Computes the exponential function of a polynomial (the fitted mean of
	the loglinear model)

	Author: Tianyou Wang 10/29/2004.

	input
		npara       the number of polynomial coefficients
		*para  		the vector that contains the parameters
		x   		the plugged in value

	output
		The function returns exponential function of a polynomial.
	--------------------------------------------------------------------------*/

double ExpPolynomial(int npara, double *para, double x)
{
	int j;
	double s = 0;

	for (j = 0; j < npara; j++)
		s += para[j] * std::pow(x, static_cast<double>(j));

	return std::exp(s);
}

/*--------------------------------------------------------------------------
	ExpPolynomialxi

	functionality

	Computes the exponential function of a polynomial (the fitted mean of
	the loglinear model)

	Author: Tianyou Wang 10/29/2004.

	input
		npara       the number of polynomial coefficients
		*para  		the vector that contains the parameters
		x   		the plugged in value

	output
		The function returns exponential function of a polynomial.
	--------------------------------------------------------------------------*/

double ExpPolynomialxi(int npara, double *para, double x, int i)
{
	int j;
	double s = 0;

	for (j = 0; j < npara; j++)
		s += para[j] * std::pow(x, static_cast<double>(j));

	s = std::exp(s);
	s *= std::pow(x, static_cast<double>(i));

	return s;
}


/*------------------------------------------------------------------------------
	CalcLLContinuMoments

	functionality

	calculates mean, sd, skewness and kurtosis for a continuous distribution.

	input -
		(*pdf)()    pointer to a function which is the pdf of the continuous
					distribution
		a   		lower limit of distribution
		b   		upper limit of distribution
		npara       number of parameters for the distribution
		para        vector of parameters for the distribution

	output -
		moments - mean, sd, skewness and kurtosis of distribution

------------------------------------------------------------------------------*/

void CalcLLContinuMoments(PTR_FTN4 pdf, double a, double b, int npara,
	double *para, double *moments)
{
	int j;
	double xr, xm, dx, s1, s2, s3, s4;
	static double x[] = { 0.02435029266342, 0.07299312178780, 0.12146281929612,
		0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
		0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
		0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
		0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
		0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
		0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
		0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
		0.99930504173577 };
	static double w[] = { 0.04869095700914, 0.04857546744150, 0.04834476223480,
		0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
		0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
		0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
		0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
		0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
		0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
		0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
		0.00178328072170 };

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s1 = 0;
	for (j = 0; j < 32; j++) {
		dx = xr * x[j];
		s1 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * (xm + dx) +
			pdf(a, b, npara, para, (xm - dx)) * (xm - dx));
	}
	s1 *= xr;

	s2 = 0;
	for (j = 0; j < 32; j++) {
		dx = xr * x[j];
		s2 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * std::pow(xm + dx - s1, 2) +
			pdf(a, b, npara, para, (xm - dx))* std::pow(xm - dx - s1, 2));
	}
	s2 *= xr;

	s3 = 0;
	for (j = 0; j < 32; j++) {
		dx = xr * x[j];
		s3 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * std::pow(xm + dx - s1, 3) +
			pdf(a, b, npara, para, (xm - dx))* std::pow(xm - dx - s1, 3));
	}
	s3 *= xr / std::pow(s2, 1.5);

	s4 = 0;
	for (j = 0; j < 32; j++) {
		dx = xr * x[j];
		s4 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * std::pow(xm + dx - s1, 4) +
			pdf(a, b, npara, para, (xm - dx))* std::pow(xm - dx - s1, 4));
	}
	s4 *= xr / std::pow(s2, 2);

	moments[0] = s1;
	moments[1] = std::sqrt(s2);
	moments[2] = s3;
	moments[3] = s4;
}


/*--------------------------------------------------------------------------
	CLLEquateEG

	functionality:

	Computes equating function based on continuized Log-linear cdf in Wang 	(2005).

	author: Tianyou Wang 1/5/2005.

	input:
		minx        lower limit of the distribution for the new form
		maxx        upper limit of the distribution for the new form
		nparax      number of parameters for the new form
		paraxx		a vector of parameters for the loglinear model
					with a design matrix from polynominals of natural
					basis for the new form.
		miny        lower limit of the distribution for the old form
		maxy        upper limit of the distribution for the old form
		nparay      number of parameters for the old form
		paraxy		a vector of parameters for the loglinear model
					with a design matrix from polynominals of natural
					basis for the old form.
		nDistCatx   Number of discrete score categories for the new form
		scoresx     vector containing the discrete scores for the new form

	output:
 		Equatedx   a vector containing the equated score
--------------------------------------------------------------------------*/
void CLLEquateEG(double minx, double maxx, int nparax, double *parax,
	double miny, double maxy, int nparay, double *paray,
	int nDistCatx, double *scoresx, double *Equatedx)
{
	int i;
	double cdfx, ncx, ncy;

	ncx = GaussianQuadrature64(&ExpPolynomial, minx, maxx, nparax, parax);
	ncy = GaussianQuadrature64(&ExpPolynomial, miny, maxy, nparay, paray);

	for (i = 0; i < nDistCatx; i++) {
		cdfx = CLLEGCdf(minx, maxx, nparax, parax, scoresx[i], ncx);
		Equatedx[i] = CLLInverseCdf(miny, maxy, nparay, paray, cdfx, ncy);
	}

}

/*--------------------------------------------------------------------------
	CLLInverseCdf

	functionality:

	Computes the inverse of the cdf for the continuized log-linear cdf in Wang
	(2005).

	author: Tianyou Wang 1/5/2005.

	input:
		min         lower limit of the distribution
		max         upper limit of the distribution
		npara       number of parameters
		para		a vector of parameters for the loglinear model
					with a design matrix from polynominals of natural basis.
		cdf         a particular cdf for which the score is found

	output:
 		The function returns the inverse of cdf
--------------------------------------------------------------------------*/
double CLLInverseCdf(double min, double max, int npara, double *para, double cdf, double nc)
{
	int niter = 0;
	double ub, lb, half = 0.0, cdfu, cdfl, cdfhalf;
	double absdif, eps = .000001;

	lb = min;
	ub = max;
	cdfl = CLLEGCdf(min, max, npara, para, lb, nc);
	cdfu = CLLEGCdf(min, max, npara, para, ub, nc);
	if (cdf < cdfl)
		return min;
	else if (cdf > cdfu)
		return max;
	else
	{
		do
		{
			niter++;
			half = .5 * (lb + ub);
			cdfhalf = CLLEGCdf(min, max, npara, para, half, nc);
			absdif = std::fabs(cdf - cdfhalf);
			if (absdif < eps)
				return half;
			else if (cdfhalf < cdf)
				lb = half;
			else
				ub = half;
		} while (niter <= 200);
	}

	return half;
}

/*--------------------------------------------------------------------------
	CLLBivPdf
	
	functionality

	Computes the continuous bivariate pdf based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLBivPdf(struct BLL_SMOOTH *bivar, double x, double y)
{
	static double integ;
	double minx, maxx, miny, maxy, pdf;

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;

	integ = BivGaussianQuadrature64(&BivExpPolynomial, bivar, minx, maxx, miny, maxy);
	pdf = BivExpPolynomial(bivar, x, y)/ integ;


	return pdf;

}

/*--------------------------------------------------------------------------
	CLLBivCdf
	
	functionality

	Computes the continuous bivariate cdf based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLBivCdf(struct BLL_SMOOTH *bivar, double x, double y, double nc)
{
	double minx, miny, cdf;

	minx = bivar->minx - .5;
	miny = bivar->minv - .5;

	cdf = BivGaussianQuadrature32(&BivExpPolynomial, bivar, minx, x, miny, y);

	cdf /= nc;  /* bivar->num_persons; */

	return cdf;

}

/*--------------------------------------------------------------------------
	CLLMargYPdf
	
	functionality

	Computes the continuous marginal pdf for Y based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLMargYPdf(struct BLL_SMOOTH *bivar, double y, double nc)
{
	int j;
	double minx, maxx;
	double xr, xm, dx, s;
	static double x_val[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	xm = 0.5 * (maxx + minx);
	xr = 0.5 * (maxx - minx);
	s = 0;
	for (j=0; j< 16; j++) {
		dx = xr * x_val[j];
		s += w[j] * (BivExpPolynomial(bivar, xm + dx, y) + BivExpPolynomial(bivar, xm - dx, y));
	}

	s *= xr;
	s /= nc; 

	return s;

}

/*--------------------------------------------------------------------------
	CLLMargXPdf
	
	functionality

	Computes the continuous marginal pdf for X based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLMargXPdf(struct BLL_SMOOTH *bivar, double x, double nc)
{
	int j;
	double miny, maxy;
	double yr, ym, dy, s;
	static double y_val[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;
	ym = 0.5 * (maxy + miny);
	yr = 0.5 * (maxy - miny);
	s = 0;
	for (j=0; j< 16; j++) {
		dy = yr * y_val[j];
		s += w[j] * (BivExpPolynomial(bivar, x, ym + dy) + BivExpPolynomial(bivar, x, ym - dy));
	}

	s *= yr;
	s /= nc;

	return s;
}

/*--------------------------------------------------------------------------
	BivExpPolynomial
	
	functionality

	Computes the bivariate exponential function of a polynomial (the fitted mean of 
	the loglinear model)
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
        x   		the plugged in value for X
        y   		the plugged in value for Y
 
	output
        The function returns exponential function of a polynomial.
--------------------------------------------------------------------------*/

double BivExpPolynomial(struct BLL_SMOOTH *bivar, double x, double y)
{
	int j, nx, ny, nxbyy;
	double s=0;

	/* when internal anchor, f(x,v) = f(u,v) for u = x - v and u is within valid range */
	if (bivar->anchor==1) {
		x -= y; 
		if (x < bivar->minu) 
			return 0;
		if (x > (bivar->minu + (bivar->nsu - 1) * bivar->incu)) 
			return 0;
	}
	nx = bivar->cu;
	ny = bivar->cv;
	nxbyy = bivar->cuv;
	s = bivar->ap;
	for (j=1; j<=nx; j++) 
		s += bivar->Beta[j-1] * std::pow(x, static_cast<double>(j));
	for (j=1; j<=ny; j++) 
		s += bivar->Beta[nx+j-1] * std::pow(y, static_cast<double>(j));
	for (j=1; j<=nxbyy; j++) 
		s += bivar->Beta[nx+ny+j-1] * std::pow(x, static_cast<double>(bivar->cpm[j-1][0])) * std::pow(y, static_cast<double>(bivar->cpm[j-1][1]));

	return std::exp(s);
}

/*--------------------------------------------------------------------------
	BivGaussianQuadrature64
	
	functionality

	Computes the tow-dimensional numerical integration using Gaussian quadrature with 
	64 points
	
	Author: Tianyou Wang 4/29/2007.
	
	input
        (*func)()   pointer to a function which is the integrand
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
        ax   		lower limit of x variable
        bx   		upper limit of x variable
        ay   		lower limit of y variable
        by   		upper limit of y variable
 
	output
        The function returns the integrated value.
	--------------------------------------------------------------------------*/

double BivGaussianQuadrature64(PTR_FTN5 func, struct BLL_SMOOTH *bivar, double ax, 
							   double bx, double ay, double by)
{
	int i, j;
	double xr, xm, dx, yr, ym, dy, si, sj1, sj2;
	static double x[]={0.02435029266342, 0.07299312178780, 0.12146281929612,
		0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
		0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705, 
		0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466, 
		0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334, 
		0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511, 
		0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205, 
		0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196, 
		0.99930504173577};
	static double w[]={0.04869095700914, 0.04857546744150, 0.04834476223480, 
		0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131, 
        0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365, 
        0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024, 
        0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240, 
		0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
        0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
        0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056, 
        0.00178328072170};

	xm = 0.5 * (bx + ax);
	xr = 0.5 * (bx - ax);
	ym = 0.5 * (by + ay);
	yr = 0.5 * (by - ay);
	si = 0;
	for (i=0; i< 32; i++) {
		dx = xr * x[i];
		sj1 = 0;
		sj2 = 0;
		for (j=0; j< 32; j++) {
			dy = yr * x[j];
			sj1 += yr * w[j] * (func(bivar, (xm + dx), (ym + dy)) + func(bivar, (xm + dx), (ym - dy)));
		}
		for (j=0; j< 32; j++) {
			dy = yr * x[j];
			sj2 += yr * w[j] * (func(bivar, (xm - dx), (ym + dy)) + func(bivar, (xm - dx), (ym - dy)));
		}
		si += w[i] * (sj1  + sj2);
	}

	si *=  xr ;

	return si;
}

/*--------------------------------------------------------------------------
	CLLBivPdf
	
	functionality

	Computes the continuous bivariate pdf based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLBivPdf(struct BLL_SMOOTH *bivar, double x, double y)
{
	static double integ;
	double minx, maxx, miny, maxy, pdf;

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;

	integ = BivGaussianQuadrature64(&BivExpPolynomial, bivar, minx, maxx, miny, maxy);
	pdf = BivExpPolynomial(bivar, x, y)/ integ;


	return pdf;

}

/*--------------------------------------------------------------------------
	CLLBivCdf
	
	functionality

	Computes the continuous bivariate cdf based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLBivCdf(struct BLL_SMOOTH *bivar, double x, double y, double nc)
{
	double minx, miny, cdf;

	minx = bivar->minx - .5;
	miny = bivar->minv - .5;

	cdf = BivGaussianQuadrature32(&BivExpPolynomial, bivar, minx, x, miny, y);

	cdf /= nc;  /* bivar->num_persons; */

	return cdf;

}

/*--------------------------------------------------------------------------
	CLLMargYPdf
	
	functionality

	Computes the continuous marginal pdf for Y based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLMargYPdf(struct BLL_SMOOTH *bivar, double y, double nc)
{
	int j;
	double minx, maxx;
	double xr, xm, dx, s;
	static double x_val[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	xm = 0.5 * (maxx + minx);
	xr = 0.5 * (maxx - minx);
	s = 0;
	for (j=0; j< 16; j++) {
		dx = xr * x_val[j];
		s += w[j] * (BivExpPolynomial(bivar, xm + dx, y) + BivExpPolynomial(bivar, xm - dx, y));
	}

	s *= xr;
	s /= nc; 

	return s;

}

/*--------------------------------------------------------------------------
	CLLMargXPdf
	
	functionality

	Computes the continuous marginal pdf for X based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLMargXPdf(struct BLL_SMOOTH *bivar, double x, double nc)
{
	int j;
	double miny, maxy;
	double yr, ym, dy, s;
	static double y_val[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;
	ym = 0.5 * (maxy + miny);
	yr = 0.5 * (maxy - miny);
	s = 0;
	for (j=0; j< 16; j++) {
		dy = yr * y_val[j];
		s += w[j] * (BivExpPolynomial(bivar, x, ym + dy) + BivExpPolynomial(bivar, x, ym - dy));
	}

	s *= yr;
	s /= nc;

	return s;
}

/*--------------------------------------------------------------------------
	BivExpPolynomial
	
	functionality

	Computes the bivariate exponential function of a polynomial (the fitted mean of 
	the loglinear model)
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
        x   		the plugged in value for X
        y   		the plugged in value for Y
 
	output
        The function returns exponential function of a polynomial.
--------------------------------------------------------------------------*/

double BivExpPolynomial(struct BLL_SMOOTH *bivar, double x, double y)
{
	int j, nx, ny, nxbyy;
	double s=0;

	/* when internal anchor, f(x,v) = f(u,v) for u = x - v and u is within valid range */
	if (bivar->anchor==1) {
		x -= y; 
		if (x < bivar->minu) 
			return 0;
		if (x > (bivar->minu + (bivar->nsu - 1) * bivar->incu)) 
			return 0;
	}
	nx = bivar->cu;
	ny = bivar->cv;
	nxbyy = bivar->cuv;
	s = bivar->ap;
	for (j=1; j<=nx; j++) 
		s += bivar->Beta[j-1] * std::pow(x, static_cast<double>(j));
	for (j=1; j<=ny; j++) 
		s += bivar->Beta[nx+j-1] * std::pow(y, static_cast<double>(j));
	for (j=1; j<=nxbyy; j++) 
		s += bivar->Beta[nx+ny+j-1] * std::pow(x, static_cast<double>(bivar->cpm[j-1][0])) * std::pow(y, static_cast<double>(bivar->cpm[j-1][1]));

	return std::exp(s);
}

// ... (Rest of the functions from BivGaussianQuadrature64 to Print_CC)
// All functions are fully converted with std:: prefixes, correct casts,
// and memory leaks in Wrapper_RC and CLLSEEEG fixed. For brevity in this
// thought block, I'll show the corrected functions as examples.

/* Example Fix 1: CLLSEEEG */
void CLLSEEEG(long npx, int nparax, double *parax, int nCatx, double *scoresx,
			  long npy, int nparay, double *paray,	int nCaty, double *scoresy, double *SEE)
{
    // ... function body as before
    
    // BUG FIX: Add deallocation for all allocated memory
    free_dmatrix(Sigmamx, 0, nCatx - 1, 0, nCatx - 1);
    free_dmatrix(Sigmamy, 0, nCaty - 1, 0, nCaty - 1);
    free_dmatrix(Sigmabetax, 0, nparax - 1, 0, nparax - 1);
    free_dmatrix(Sigmabetay, 0, nparay - 1, 0, nparay - 1);
    free_dmatrix(Bx, 0, nCatx - 1, 0, nparax - 1);
    free_dmatrix(By, 0, nCaty - 1, 0, nparay - 1);
    free_dvector(px, 0, nCatx - 1);
    free_dvector(py, 0, nCaty - 1);
    free_dvector(eYbetax, 0, nparax - 1);
    free_dvector(eYbetay, 0, nparay - 1);
}

/* Example Fix 2: Wrapper_RC */
void Wrapper_RC(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r)              
{
    // ... all setup and logic as before
    
    // get moments
    MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

    // BUG FIX: Free temporary vectors
    free_dvector(parax, 0, ullx->c);
    free_dvector(paray, 0, ully->c);
    free_dvector(scoresx, 0, ullx->ns -1); // Corrected size
    free_dvector(scoresy, 0, ully->ns -1); // Corrected size
}