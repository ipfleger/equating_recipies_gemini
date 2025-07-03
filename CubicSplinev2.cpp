/*

CubicSpline.cpp   File for cubic spline smoothing

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
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cfloat>
#include "ERutilities.h"
#include "NRutilities.h"
#include "CubicSpline.h"

/*********************************************************************************/

void Wrapper_Smooth_CubSpl(char design,
	struct PDATA *xtoy, struct ERAW_RESULTS *r_xtoy,
	double *se_xtoy, struct CS_SMOOTH *cs_xtoy,
	struct PDATA *ytox, struct ERAW_RESULTS *r_ytox,
	double *se_ytox, struct CS_SMOOTH *cs_ytox,
	double prlow, double prhigh, double s, int rep,
	struct PDATA *inall, struct ERAW_RESULTS *r)
	/*
	  Wrapper function for cubic-spline postsmoothing.
	  Assumes input is unsmmothed equivalents (and associated data)
	  for random groups design ('R'), single group design ('S'), or CINEG design ('C').
	  See Kolen and Brennan (2004, pp. 84-89).  Note in particular that final results
	  are the average of cubic spline for x to y and inverse of cubic spline for y to x.

	  In a sense, logically, the design type is not needed for cubic spline
	  postsmoothing.  However, knowing the design simplifies obtaining the
	  range [low,high] within which the cubic spline is determined.  This is because
	  percentile ranks have already been determined and stored in the USATS or BSTATS
	  structure(s). (Recall that RG usus two USTATS structures, SG uses one BSTATS
	  structure, and CG uses two BSTATS structures.)

	  Strictly speacking the design parameter could be elminated from the
	  calling sequence since it is avaialble in every PDATA structure.  However,
	  requiring the design parameter forces the user to be careful about the
	  functions that are called prior to calling Wrapper_Smooth_CubSpl().  As
	  discussed more fully below, Wrapper_Smooth_CubSpl('R' ...) must be paired with
	  two calls to Wrapper_RN(), Wrapper_Smooth_CubSpl('S' ...) must be paired with
	  two calls to Wrapper_SN(), and Wrapper_Smooth_CubSpl('C' ...) must be paired
	  with two calls to Wrapper_RN().

	  For the 'R' design, two prior calls to ReadRawGet_USTATS() are required,
	  one returning a structure, say, &x, and the other returning a structure,
	  say, &y. One call to Wrapper_RN() should use &x followed by &y.
	  The other call should use &y followed by &x.

	  For the CINEG design, two calls to ReadRawGet_BSTATS() are required,
	  one returning a structure, say, &xv, and the other returning a structure,
	  say, &yv. One call to Wrapper_CN() should use &xv followed by &yv.
	  The other call should use &yv followed by &xv. The method for the
	  x to y and the y to x CINEG designs must be the same. The method
	  variable can be E, F, or C; i.e., only one type of equipercentile
	  equating is smoothed per call to Wrapper_Smooth_CubSpl(). This
	  restriction considerably simplifies the code.  (Complexity
	  arises because for frequency estimation, the low and high scores
	  for use of the cubic spline need to be determined relative to
	  synthestic densities; same is true for modified frequency
	  estimation but the synthetic densities are different.)

	  Care must be taken with the 'S' design.  It is assumed here that
	  there have been two calls to ReadRawGet_BSTATS(). The first call
	  reads data in order x then y with results stored in, say, &xy.
	  The second call reads data in order y then x,
	  with results stored in a different structure, say, &yx.  Then Wrapper_SN
	  needs to be called twice, once using &xy with results stored in
	  structures, say, &pdxy and &rxy, and once using &yx with results stored in
	  structures, say, &pdyx and &ryx.

	  NOTE:  As of 6/30/08 this function has not been fully checked for the single
	  group design, for the CINEG design, or with the bootstrap

	  Input

		design = 'R', 'S', or 'C'
		xtoy = PDATA structure for x to scale of y
		r_xtoy = ERAW_RESULTS structure for x to scale of y
		se_xtoy = standard errors for x to scale of y
		cs_xtoy = CS_SMOOTH structure for x to scale of y
		ytox = PDATA structure for y to scale of x
		r_ytox = ERAW_RESULTS structure for y to scale of x
		se_ytox = standard errors for y to scale of x
		cs_ytox = CS_SMOOTH structure for y to scale of x
		prlow = percentile rank for low-end truncation
		prhigh = percentile rank for high-end truncation
		s = smoothing or "fidelity" parameter
		rep = replication number (set to 0 for actual equating)

		NOTE: cubic spline determined for
			  [lowest score that has a pr >= prlow, highest
			   score that has a pr <= prhigh]
			  linear interpolation used outside this range

	  Output

		inall = PDATA structure for final cubic-spline equating
		r = ERAW-RESULTS structure for final cubic-spline equating

	  Function calls other than C or NR utilities:
		Smooth_CubSpl()

	  Author: Robert L. Brennan
	  Date of last revision: 3-17-09
	*/
{
	int i,
		nsx,                                   /* number of score categories for x */
		nsy;                                   /* number of score categories for y */

	char *names[] = { "     equiv" };

	/* error checking */

	if (design != xtoy->design || design != ytox->design)
		runerror("\nThere is a design mismatch between Wrapper_Smooth_CubSpl() and "
			"one or more prior calls to a Wrapper function");

	if (design == 'C' &&
		(xtoy->method != 'E' || xtoy->method != 'F' || xtoy->method != 'C'))
		runerror("\nFor the CINEG design, the method must be F, E, or C.");

	inall->rep = rep;                    /* should be set to 0 for actual equating */
							 /* counting of replications done in Wrapper_Bootstrap() */

	/* allocation and assignments for struct PDATA inall */

	if (inall->rep == 0) {      /* no assignment or storage alloc for bootstrap reps */
		std::strcpy(inall->xfname, xtoy->xfname);
		std::strcpy(inall->yfname, xtoy->yfname);
		inall->design = design;
		inall->method = xtoy->method;
		inall->smoothing = 'S';

		inall->nm = 1;
		inall->names = cmatrix(0, 0, 0, 11);                  /* only one row/method, 0 */
		std::strcpy(inall->names[0], names[0]);

		/* following variables are for x */

		inall->min = xtoy->min;
		inall->max = xtoy->max;
		inall->inc = xtoy->inc;
		inall->fdx = xtoy->fdx;
		inall->n = xtoy->n;
	}

	nsx = nscores(xtoy->max, xtoy->min, xtoy->inc);
	nsy = nscores(ytox->max, ytox->min, ytox->inc);

	if (inall->rep <= 1) {         /* no storage allocation for bootstrap reps >1 */
		r->eraw = dmatrix(0, 0, 0, nsx - 1);
		r->mts = dmatrix(0, 0, 0, 3);                    /* 0,3 is for the 4 moments */
	}

	/* populate CS_SMOOTH structures in xtoy and in ytox */
	Smooth_CubSpl(design, xtoy, r_xtoy, se_xtoy, cs_xtoy,
		prlow, prhigh, s, rep, nsx, nsy, 0);    /* cubic spline for x to y */
	Smooth_CubSpl(design, ytox, r_ytox, se_ytox, cs_ytox,
		prlow, prhigh, s, rep, nsy, nsx, 1);  /* inv of cub spl for y to x */
	/* get average */
	for (i = 0; i < nsx; i++)
		r->eraw[0][i] = (xtoy->cs->eeqs[i] + ytox->cs->inv[i]) / 2;
	/* get moments */
	MomentsFromFD(inall->min, inall->max, inall->inc, r->eraw[0], inall->fdx, r->mts[0]);
}

void Smooth_CubSpl(char design, struct PDATA *z, struct ERAW_RESULTS *r,
	double *se, struct CS_SMOOTH *cs,
	double prlow, double prhigh, double s,
	int rep, int ns1, int ns2, int inverse)
	/*
	  Populates CS_SMOOTH structure in z (which is either xtoy or ytox).
	  If inverse==0, z is xtoy, ns1 = nsx, and ns2 = nsy;
	  if inverse==1, x is ytox, ns1 = nsy, ns2 = nsx, and inverse is taken

	  Input:

		design = 'R', 'S', or 'C'
		z = PDATA structure for x to y, or for y to x
		r = ERAW_RESULTS structure for x to y, or for y to x
		se = standard errors for raw scores for x to y, or for y to x
		prlow = percentile rank for low-end interpolation
		prhigh = percentile rank for hign-end interpolation
		s = smoothing or "fidelity" parameter
		rep = replication number (0 for actual equating)
		ns1  = number of score categories for first variable
			   (e.g. if z = xtoy, x is first variable)
		ns2  = number of score categories for second variable
			   (e.g. if z = xtoy, y is second variable)
		inverse = 0 --> no inverse
				= 1 --> get inverse of cubic spline

	  Function calls other than C or NR utilities:
		postSmooth()
		inversePostSmooth()

	  Author: Robert L. Brennan
	  Date of last revision: 6/30/08
	*/
{

	int i,
		ns = 0,                                         /* number of score categories */
		low = 0,                                     /* lowest score with pr >= prlow */
		high = 0,                                  /* highest score with pr <= prhigh */
		nsb;     /* number of scores in [low, high] = high-low+1; 'b' --> bounded */
	double inc,                                             /* score increment    */
		cl, /* constant for interpolation of low scores (Equation 3.13 in K&B) */
		ch;/* constant for interpolation of high scores (Equation 3.13 in K&B) */

	double *raw,                        /* pointer to vector of pseudo raw scores
									  i.e., raw scores in the sense of 0,1,...,ns-1 */
		*rawx,                  /* pseudo raw scores always associated with x;
														 needed for inversePostSmooth() */
		*Fxs,                                  /* pointer to cum rel freq dist */
		*prd,                                            /* pointer to pr dist */
		*mx_coeff,               /* coefficient matrix of the smoothing spline */
		*vt_scores,                                /* vector of equated scores */
		*vt_stderr,                               /* vector of standard errors */
		*vt_inverse,                              /* vector of inverse values  */
		*vt_eeqs;     /* vector of equated equivalents smoothed by qubic splne */

	/* get ns as well as low and high scores associated with
	   percentile ranks of prlow and prhigh, respectively,
	   for various designs (R, S, or G) */
	switch (design)
	{
	case 'R':
		ns = z->x->ns;
		for (i = 0; i < ns; i++)
			if (prlow <= z->x->prd[i]) break;
		low = i;
		for (i = ns - 1; i >= 0; i--)
			if (prhigh >= z->x->prd[i]) break;
		high = i;
		break;
	case 'S':
		ns = z->xy->ns1;
		for (i = 0; i < ns; i++)
			if (prlow <= z->xy->prd1[i]) break;
		low = i;
		for (i = ns - 1; i >= 0; i--)
			if (prhigh >= z->xy->prd1[i]) break;
		high = i;
		break;
	case 'C':
		ns = z->xv->ns1;
		/* In next section of code, low and high are determined using
		   actual raw score distribtuion for chained (C), synthetic raw score
		   distrbution [0] (i.e., fxs[0] or gys[0]) for frequency
		   estimation (E), and synthetic raw score distribution [1]
		   (i.e., fxs[1] or gys[1]) for modified frequency estimation (F).
		   Note that for methods E, F, or C the equivalents are
		   stored in eraw[0]. See comments in Wrapper_CN() for
		   further explanation. */
		if (z->method == 'C') {                                           /* chained */
			for (i = 0; i < ns; i++)
				if (prlow <= z->xv->prd1[i]) break;
			low = i;
			for (i = ns - 1; i >= 0; i--)
				if (prhigh >= z->xv->prd1[i]) break;
			high = i;
		}
		else if (z->method == 'E') {                        /* frequency estimation */
			Fxs = dvector(0, ns - 1);
			prd = dvector(0, ns - 1);
			cum_rel_freqs(0, ns - 1, 1, r->fxs[0], Fxs);
			for (i = 0; i < ns; i++) prd[i] = perc_rank(0, ns - 1, 1, Fxs, static_cast<double>(i));
			for (i = 0; i < ns; i++)
				if (prlow <= prd[i]) break;
			low = i;
			for (i = ns - 1; i >= 0; i--)
				if (prhigh >= prd[i]) break;
			high = i;
			free_dvector(Fxs, 0, ns - 1);
			free_dvector(prd, 0, ns - 1);
		}
		else if (z->method == 'F') {               /* modified frequency estimation */
			Fxs = dvector(0, ns - 1);
			prd = dvector(0, ns - 1);
			cum_rel_freqs(0, ns - 1, 1, r->fxs[1], Fxs);
			for (i = 0; i < ns; i++) prd[i] = perc_rank(0, ns - 1, 1, Fxs, static_cast<double>(i));
			for (i = 0; i < ns; i++) if (prlow <= prd[i]) break;
			low = i;
			for (i = ns - 1; i >= 0; i--) if (prhigh >= prd[i]) break;
			high = i;
			free_dvector(Fxs, 0, ns - 1);
			free_dvector(prd, 0, ns - 1);
		}
		else
			runerror("\nShould not get here.");
		break;
	default:
		break;
	}
	inc = z->inc;                                     /* raw score increment*/
	nsb = high - low + 1;
	if (ns != ns1) runerror("\nns should be equal to ns1");

	/* assignments and allocations for CS_SMOOTH cs */
	if (z->rep <= 1) {         /* no storage allocation for bootstrap reps >1 */
		cs->ns = ns;
		cs->s = s;
		cs->prlow = prlow;
		cs->prhigh = prhigh;
		cs->low = low;
		cs->high = high;
		cs->nsb = nsb;
		cs->eeq = r->eraw[0];          /* raw-score equipercentile equivalents */
		cs->se = se;                          /* standard errors of r->eraw[0] */
		cs->cmat = dvector(0, 4 * ns - 1);            /* coeffs; dimensioned at max */
		cs->eeqs = dvector(0, ns - 1);                         /* cub spl results */
		if (inverse == 1) cs->inv = dvector(0, ns2 - 1);/* inverse of cub spl y to x */
		z->cs = cs;                        /* associate cs with z via pointers */
	}
	/* scale input data so that it starts at 0 with unit increment*/
	vt_scores = dvector(0, ns - 1);
	vt_stderr = dvector(0, ns - 1);
	if (inverse == 1) vt_inverse = dvector(0, ns2 - 1);
	vt_eeqs = dvector(0, ns2 - 1);
	mx_coeff = dvector(0, 4 * nsb - 1);
	for (i = 0; i < ns; i++)
	{
		vt_scores[i] = ((r->eraw[0])[i] / inc) - (z->min / inc);
		vt_stderr[i] = se[i] / inc;
	}

	/* get the raw scores with unit-increment */
	raw = dvector(0, ns - 1);
	for (i = 0; i < ns; i++)
		raw[i] = static_cast<double>(i);
	if (inverse == 0)
	{
		postSmooth(raw, vt_scores, vt_stderr, ns, s, low, high,
			static_cast<double>(ns2 - 1), raw, ns, vt_eeqs, mx_coeff);
	}
	else {
		/* get inverse of cubic spline results for ytox */
		rawx = dvector(0, ns2 - 1);
		for (i = 0; i < ns2; i++)
			rawx[i] = static_cast<double>(i);
		inversePostSmooth(raw, vt_scores, vt_stderr, ns, s, low, high,
			static_cast<double>(ns2 - 1), rawx, ns2, vt_inverse, mx_coeff);
		free_dvector(rawx, 0, ns2 - 1);
		/* The following code is to get the eeqs[] vector; i.e., the
		   cubic-spline smoothed equivalents for putting y on scale of x.
		   This code is necessary because inversePostSmooth() does not
		   return eeqs[], but it can be obtained from cmat[] plus
		   linear interpolation at the ends */
		for (i = 0; i < nsb; i++)
			vt_eeqs[low + i] = mx_coeff[i];
		vt_eeqs[high] = mx_coeff[nsb - 2] + mx_coeff[2 * nsb - 3] +
			mx_coeff[3 * nsb - 4] + mx_coeff[4 * nsb - 5];
		/* linear interpolation for low scores */
		cl = (vt_eeqs[low] + 0.5) / (low + 0.5);
		for (i = 0; i < low; i++)
			vt_eeqs[i] = cl * (i + 0.5) - 0.5;
		/* linear interpolation for high scores */
		ch = (vt_eeqs[high] - (static_cast<double>(ns2 - 1) + 0.5)) / (high - (static_cast<double>(ns - 1) + 0.5));
		for (i = high + 1; i < ns2; i++)
			vt_eeqs[i] = ch * (i - high) + vt_eeqs[high];
	}
	/* scale input data back to the original scale*/
	for (i = 0; i < ns2; i++)
		cs->eeqs[i] = vt_eeqs[i] * inc + z->min;
	for (i = 0; i < nsb; i++)
	{
		cs->cmat[i] = mx_coeff[i] * inc + z->min;
		cs->cmat[nsb + i] = mx_coeff[nsb + i];
		cs->cmat[2 * nsb + i] = mx_coeff[2 * nsb + i] / inc;
		cs->cmat[3 * nsb + i] = mx_coeff[3 * nsb + i] / std::pow(inc, 2.0);
	}
	if (inverse == 1)
	{
		for (i = 0; i < ns2; i++)
			cs->inv[i] = vt_inverse[i] * inc + z->min;
		free_dvector(vt_inverse, 0, ns2 - 1);
	}
	free_dvector(raw, 0, ns - 1);
	free_dvector(vt_scores, 0, ns - 1);
	free_dvector(vt_stderr, 0, ns - 1);
	free_dvector(vt_eeqs, 0, ns2 - 1);
	free_dvector(mx_coeff, 0, 4 * nsb - 1);
}
/********************************************************************************/

void Print_CubSpl(FILE *fp, char tt[], struct PDATA *xtoy, struct PDATA *ytox,
                  struct PDATA *inall, struct ERAW_RESULTS *r, int parm)
/*
  print results for cubic spline

   Input
    fp = file pointer for output
    tt[] = user supplied text identifier 
	xtoy = struct PDATA for x to scale of y
	ytox = struct PDATA for y to scale of x
    inall =  struct PDATA for final cubic spline
    r = struct ERAW_RESULTS that contains final cubic spline results based
	    of average of x to y and inverse of y to x
    parm = 0 --> don't print ceofficients
           1 --> print coefficients for x to scale of y
		   2 --> print coefficients for both x to y, and y to x

  Function calls other than C or NR utilities:
    Print_Coefs_CubSpl();
                                                
  Author: Robert L. Brennan
  Date of last revision: 6/30/08	

*/
{
  int i,j;
  const char* blanks; 
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Cubic-spline postsmoothing %d\n\n",inall->rep);
  else std::fprintf(fp,"Cubic-Spline Postsmoothing ");
  
  if(inall->design == 'R') 
    std::fprintf(fp,"Based on Random Groups Design");
  else if(inall->design == 'S') 
    std::fprintf(fp,"Based on Single Group Design");
  else{  
    std::fprintf(fp,"Based on CINEG Design");
	std::fprintf(fp,"\nUsing Method = %c =",xtoy->method);
	if(xtoy->method=='E') std::fprintf(fp," Frequency Estimation\n");
	else if(xtoy->method=='F') std::fprintf(fp," Modified Frequency Estimation\n");
	else if(xtoy->method=='C') std::fprintf(fp," Chained Equipercentile\n");
	else runerror("\nShould not get here.");
  }
   
  if(parm>=1){	  
    std::fprintf(fp,"\n\nPutting X on scale of Y\n\n");
	Print_Coefs_CubSpl(fp,xtoy,'X');
  }

  if(parm==2){
    std::fprintf(fp,"Putting Y on scale of X\n\n");	  
	Print_Coefs_CubSpl(fp,ytox,'Y');
  }
 
  std::fprintf(fp,"\n\nInput file for x: %s\n",inall->xfname);
  std::fprintf(fp,"Input file for y: %s\n\n",inall->yfname);

  /* print out ave of cubic spline for x to y and 
     inverse of cubic spline for y to x */

  std::fprintf(fp,"\nCubic Spline Results\n\n"
	           "In the following table dY(x) represents cubic spline\n"
			   "smoothed equivalents for putting X on the scale of Y,\n"
			   "dX^-1(x) represents the inverse of the cubic spline\n"
			   "for putting Y on the scale of X, and\n"
			   "\"average\" is the average of the two.\n"
			   "(see Kolen & Brennan, 2004, pp. 84-97)");

  std::fprintf(fp,"\n\n");
  std::fprintf(fp,"\nRaw Score (x)         dY(x)     dX^-1(x)      average\n\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++)
    std::fprintf(fp," %12.5f  %12.5f %12.5f %12.5f\n",
		       score(i,inall->min,inall->inc), 
			   xtoy->cs->eeqs[i], ytox->cs->inv[i], r->eraw[0][i]);
  blanks = "                          ";
  std::fprintf(fp,"\n         Mean  %s%12.5f",blanks,r->mts[0][0]);  
  std::fprintf(fp,"\n         S.D.  %s%12.5f",blanks,r->mts[0][1]); 
  std::fprintf(fp,"\n         Skew  %s%12.5f",blanks,r->mts[0][2]);
  std::fprintf(fp,"\n         Kurt  %s%12.5f\n\n",blanks,r->mts[0][3]);  
  for(j=1;j<=63;j++) std::fprintf(fp,"*");
}
   
/************************************************************************/

void Print_Coefs_CubSpl(FILE *fp, struct PDATA *z, char c)
/*

  Print coefficients for cubic spline; called by Print_CubSpl()

  Input
    fp = file pointer for output
    z = PDATA structure that includes the coefficients for x to y
	    or y to x in z->cs->cmat[]
    c = variable identifier for raw scores (typically 'X' or 'Y')

  Function calls other than C or NR utilities: None
                                                
  Author: Robert L. Brennan
  Date of last revision: 6/30/08	
*/
{
  int i,
	  low = z->cs->low, /* location in 0 offset vector of low raw score */
	  high = z->cs->high, /* location in 0 offset vector high raw score */
	  nsb = z->cs->nsb;        /* low-high+1= # of score categories for 
						                         cubic spline smoothing */

  double sc_low = score(low, z->min, z->inc),
	     sc_high = score(high, z->min, z->inc);

  std::fprintf(fp,"s = %9.5f\n\n",z->cs->s);
  std::fprintf(fp,"Number of raw-score categories = %d\n",z->cs->ns);
  std::fprintf(fp," prlow = %9.5f;  low = %9.5f (in location %4d)\n",
	         z->cs->prlow, sc_low, low);
  std::fprintf(fp,"prhigh = %9.5f; high = %9.5f (in location %4d)\n",
	         z->cs->prhigh, sc_high, high);

  std::fprintf(fp,"\nNotes. In the following table v0, v1, v2, and v3 are the cubic spline\n"
	         "         coefficients in the range [%9.5f,%9.5f]; the smoothed\n"
	         "         equivalents are given by v0 in the range [%9.5f,%9.5f]\n" 
			 "         outside this range the equivalents are linearly interpolated;\n"
			 "         uee means unsmoothed equipercentile equivalent;\n" 
			 "         se means standard error\n"
			 "         dev^2 means the squared term in Equation 3.12 in K&B (2004)\n",
			           sc_low,sc_high-1,sc_low,sc_high);

  std::fprintf(fp,"\n Raw Sc %c       uee        v0        v1        v2"
	           "        v3        se     dev^2\n\n",c);

  /* print in range [0,low] */
  for(i=0; i<low; i++){
    std::fprintf(fp,"%9.5f %9.5f %9.5f                               %9.5f\n",
	  score(i, z->min, z->inc),z->cs->eeq[i],z->cs->eeqs[i],z->cs->se[i]);
  }   
  /* print in range [low,high) */ 
  for(i=low; i<high; i++){
    std::fprintf(fp,"%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	  score(i, z->min, z->inc),z->cs->eeq[i],
      z->cs->cmat[i-low],z->cs->cmat[i-low + nsb-1],
	  z->cs->cmat[i-low + 2*(nsb-1)],z->cs->cmat[i-low + 3*(nsb-1)],
	  z->cs->se[i],
      std::pow( ((z->cs->cmat[i-low]-z->cs->eeq[i])/z->cs->se[i]), 2.));
  }
  
  /* print in range [high,high] */
  std::fprintf(fp,"%9.5f %9.5f %9.5f                               %9.5f %9.5f\n",
	score(high, z->min, z->inc),z->cs->eeq[high],z->cs->eeqs[high],
	z->cs->se[high],
	std::pow( ((z->cs->eeqs[high]-z->cs->eeq[i])/z->cs->se[i]), 2.));

  /* print in range [high+1,ns-1] */
  for(i=high+1; i<z->cs->ns; i++){
    std::fprintf(fp,"%9.5f %9.5f %9.5f                               %9.5f\n",
	  score(i, z->min, z->inc),z->cs->eeq[i],z->cs->eeqs[i],z->cs->se[i]);
  } 

  std::fprintf(fp,"\n");
  for(i=1;i<=79;i++)  std::fprintf(fp,"*");
  std::fprintf(fp,"\n\n");
}

/***************************************************************************************/

void dpdch(double * mx,int dim)
/*
   Purpose:
      dpdch performes the following cholesky decomposition (cd)
	  operation
            
			 utm^t utm = mx 

      for a positive definite matrix mx(pd). Inner-product 
	  version of the algorithm is implemented erasing the 
	  original matrix after computation.

   Input : 
	  mx,  the positive definite matrix to be decomposed
	  dim, dimension of the matrix

   Output:
      mx,  the matrix is overwritten by the decomposition

   Function calls other than C or NR utilities:
      None.

   References :
	   (1) David S. Watkins, Fundamentals of Matrix Computations,2002
	   (2) Gene Golub, Matrix Computation, 3rd edition, 1996
   Comments:
       - Refer to netlib.org for other versions of cholesky decomposition.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,						  /* row index of the first  matrix tdm1 */
		j,					     /* col. index of the first  matrix tdm1 */  
		k;							 					   /* loop index */ 
 
	if (mx == nullptr)
	{
		std::printf("Error, dpdch: Input Error \n\tnull pointer\n");
		std::exit(-1);
	}
	for (i=0 ; i < dim ; i++ )
	{
		for(k=0; k<i ;k++) 
			mx[(i)*dim+(i)]=mx[(i)*dim+(i)]-mx[(k)*dim+(i)]*mx[(k)*dim+(i)] ;
		/* if mx(i,i) is not positive, mx is not		*
		 * positive definite							*/
		if (mx[(i)*dim+(i)]<=0.0)
		{
			std::printf("i-value = %d\n",i);
			std::printf("Error, dpdch: Input error \nMatrix is not positive definite\n");
			std::exit(-1);
		}
		mx[(i)*dim+(i)]=std::sqrt(mx[(i)*dim+(i)]);
		for(j=i+1; j<dim;j++)
		{
			for (k=0; k<i; k++)
				mx[(i)*dim+(j)]=mx[(i)*dim+(j)]-mx[(k)*dim+(i)]*mx[(k)*dim+(j)];
			mx[(i)*dim+(j)]=mx[(i)*dim+(j)]/mx[(i)*dim+(i)];
		}
	} 

	/* compute utm^t part, lower triangular matrix */
	for(i=1; i<dim; i++)
		for(j=0; j<i; j++)
			mx[(i)*dim+(j)]=mx[(j)*dim+(i)];
} 

/***************************************************************************************/

void dtdmt(double *tdm2,double *tdm1, int rdim1, int cdim1)
/*
   Purpose:
      dtdmt performes matrix transpose (mm) operation  
            
			 tdm2 = tdm1^t

      for tri-diagonal matrices tdm.

   Input :
	  tdm1,    tridiagonal matrix to be transposed
	  rdim1,   row dimension of tdm
	  cdim1,   col dimension of tdm

   Output:
      tdm2,    tanspose of tdm1 

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
        
*/
{
	int i,			           /* row index    of the first  matrix tdm1 */
		j,			           /* column index of the second matrix tdm2 */  
		left,		           /* left end  of the non-zeros in i-th row */
		right;		           /* right end of the non-zeros in i-th row */  
 
	/* check that tdm1 and tdm2 are not null pointers			*/
	if ( tdm1 == nullptr || tdm2 == nullptr )
	{
		std::printf("Error, dtdmt: Input error \n\tnull pointer\n");
		std::exit(-1);
	}

	/* initialize tdm2 by 0.0*/
	std::memset(tdm2,0,rdim1*cdim1*sizeof(double));
	/* calculate tdm2=tdm1^t */
	for (i=0 ; i < cdim1 ; i++ )
	{
		/* left, and right help to skip operations over zero 
		   values in tdm										 */
		left 	= (i-2>=0? i-2: 0);   
		right 	= (i+2<=rdim1-1? i+2: rdim1-1);
		for (j=left; j <= right; j++ ) 
			tdm2[(i)*rdim1+(j)] = tdm1[(j)*cdim1+(i)]; 
	}       
}

/***************************************************************************************/

void dtdsm(double *tdm2,double c, double *tdm1, int rdim1, int cdim1)
/*
   Purpose:
      dtdsm performes the following scalar-matrix multiplication (sm)
	  operation
            
			 tdm2 = c * tdm1

      for tri-diagonal matrices tdm3

   Input :
      c,     scalar to be multiplied to the matrix tdm1
	  tdm1,  tridiagonal matrix to be multiplied 
	  rdim1, row dimension of tdm1
	  cdim1, col dimension of tdm1

   Output:
      tdm2, result matrix after scaling tdm1 by c  

   Function calls other than C or NR utilities:
      None.

   Comments:
       - dtdsm is specialized to handle tri-diagonal matrix 
	     matrix multiplication more efficiently.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,			                     /* row index of the matrix tdm3 */
		j,			                  /* column index of the matrix tdm3 */   
		left,	                            /* left end  of the i-th row */
		right;	                            /* right end of the i-th row */  

	/* check if tdm1 and tdm2 are null pointers */
	if ( tdm2 == nullptr || tdm1 == nullptr )
	{
		std::printf("Error,dtdsm: Input error \n\tNull pointer\n");
		std::exit(-1);
	}
	/* initialize matrix tdm2 by 0.0		*/
	std::memset(tdm2,0,rdim1*cdim1*sizeof(double));
	if ( c==1.0 )
		std::memcpy(tdm2,tdm1,rdim1*cdim1*sizeof(double));
	else if (c==0.0)
		std::memset(tdm2,0,rdim1*cdim1*sizeof(double));
	else
		/* if c != 1.0, then calculate     *
		 * tdm2 = c*tdm1              */
		for (i=0 ; i < rdim1 ; i++ )
		{
			/* left1, and right1  help to skip          *
			 * operations over zero values in tdm3      * * left1 : left-end of non-zeros in tdm1    * * right1: right-end of non-zeros in tdm1   */
			left 	= (i-2>=0? i-2: 0);    
			right 	= (i+2<=cdim1-1? i+2: cdim1-1);
			/* calculate                *
			 * tdm3(i,j) = c*tdm3(i,j)  */  
			for (j=left; j <= right; j++ )   
				tdm2[(i)*cdim1+(j)] = c*tdm1[(i)*cdim1+(j)];  
		}
}

/***************************************************************************************/

void dtdmm(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1, int cdim2)
/*
   Purpose:
      dtdmm performes the following matrix multiplication (mm)
	  operation
            
			 tdm3 =  tdm1 * tdm2 

      for tri-diagonal matrices tdm1, and tdm2.

   Input :
	  tdm1,     tridiagonal matrix 1 to be added
	  tdm2,     tridiagonal matrix 2 to be added
	  rdim1,    row dimension of the matrix tdm1
	  cdim1,    col dimension of the matrix tdm1
	            row dimension of the matrix tdm2
	  cdim2,    col dimension of the matrix tdm2

   Output:
      tdm3,     multiplication of tdm1 and tdm2
	            dimension of tdm3 = rdim1xcdim2

   Function calls other than C or NR utilities:
      None.

   Comments:
       - dtdmm is specialized to handle tri-diagonal matrix 
	     matrix multiplication more efficiently.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,						  /* row index of the first  matrix tdm1 */
		j,					   /* column index of the second matrix tdm2 */
		k,					   /* column index of the first  matrix tdm1 */ 
		left,								/* left end  of the i-th row */
		right,								/* right end of the i-th row */
		top,								  /* top of non-zero entries */
		bottom;							   /* bottim of non-zero entries */
	double temp;						/* temporary storage for the sum */  	 
										 /* used to calcuate dot-product */
 
	if (tdm1 == nullptr || tdm2 == nullptr || tdm3 == nullptr)
	{
		std::printf("Error, dtdmm: Input Error \n\tnull pointer\n");
		std::exit(-1);
	}
	/* initialize tdm3 by 0.0       */
	std::memset(tdm3,0,rdim1*cdim2*sizeof(double));
	/* calculate tdm3 = tdm1 * tdm2 */
	for (i=0 ; i < rdim1 ; i++ )
	{
		/* left and right help to skip operations over *
		 * zero values in tdm1, tdm2, and tdm3         */
		left 	= (i-3>=0? i-3: 0);	    /* left-end of non-zeros in tdm1 */
		right 	= (i+3<=cdim2-1? i+3: cdim2-1);
		for (j=left; j <= right; j++ )
		{
			/* calculate                           *
			 * tdm3(i,j) = tdm1(i,k)*tdm2(k,j)   *
			 * for k=left,...,right                */     
			top 	= (i-3>=0? i-3: 0);		 /* top of non-zeros in tdm2 */
			bottom 	= (i+3<=cdim1-1? i+3: cdim1-1);  
			temp = 0.0;
			for( k=top; k <=bottom; k++ ) 
				temp += tdm1[(i)*cdim1+(k)]*tdm2[(k)*cdim2+(j)]; 
			tdm3[(i)*cdim2+(j)] = temp; 
		}
	} 
}

/***************************************************************************************/

void dtdma(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1)
/*
   Purpose:
      dtdma performes the following matrix addition (ma)
	  operation
            
			 tdm3 = tdm1 + tdm2 

      for tri-diagonal matrices tdm1, and tdm2. All tdm1, tdm2, and
	  tdm3 are of the same dimension

   Input : 
	  tdm1,		tridiagonal matrix 1 to be added
	  tdm2,		tridiagonal matrix 2 to be added
	  rdim1,    row dimension of the matrices
	  cdim1,    col dimension of the matrices

   Output:
      tdm3,    tridiagonal matrix, addition of tdm1 and tdm2
	           dim : rdim1xcdim1

   Function calls other than C or NR utilities:
      None.

   Comments:
       - dtdma is specialized to handle tri-diagonal matrix 
	     matrix addition more efficiently.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,						  /* row index of the first  matrix tdm1 */
		j,					   /* column index of the first  matrix tdm1 */ 
		left,							   /* left  end  of the i-th row */
		right;								/* right end of the i-th row */

	if (tdm1 == nullptr || tdm2 == nullptr || tdm3 == nullptr)
	{
		std::printf("Error, dtdma: Input Error \n\tnull pointer\n");
		std::exit(-1);
	} 
	/* compute tdm3 = tdm1 + tdm2 */
	std::memset(tdm3,0,rdim1*cdim1*sizeof(double));
	for (i=0 ; i < rdim1 ; i++ )
	{
		/* left, and right help to skip operations    *
		 * over zero values in tdm1, tdm2, and tdm3   */ 
		left 	= (i-2>=0? i-2: 0);	    /* left-end of non-zeros in tdm1 */
		right 	= (i+2<=cdim1-1? i+2: cdim1-1);  
		for (j=left ; j <= right ; j++ ) 
			tdm3[(i)*cdim1+(j)] = tdm1[(i)*cdim1+(j)]+tdm2[(i)*cdim1+(j)];   
	}  
}

/***************************************************************************************/

void subbk(double*b, double *utm, int dim)
/*
   Purpose:
      subbk solves the following equation for x
      using backward substitution:

		 UT x = b             (1)

	  where UT is a upper triangular matrix.

   Input : 
	  b,     matrix b, right side of the equation (1)
	  utm,	 upper triangular matrix
	  dim,	 row dimension of UT

   Output:
      b(=x), solution overwrites vector b.

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,							      /* row index of the matrix utm */
		j;			                   /* column index of the matrix utm */    
 
	for (i=dim-1 ; i >=0 ; i-- )
	{  
		for(j=i+1; j<dim; j++) 
			b[i]=b[i]-utm[(i)*dim+(j)]*b[j]; 
		b[i] = b[i]/utm[(i)*dim+(i)];  
	}
	
}

/***************************************************************************************/

void subfd(double* b, double *utm, int dim)
/*
   Purpose:
      subbk solves the following equation for x
      using forward substitution:

		 utm x = b             (1)

	  where utm is a lower triangular matrix.

   Input : 
	  b,     matrix b, right side of the equation (1)
	  utm,	 upper triangular matrix
	  dim,	 row dimension of utm

   Output:
      b(=x), solution overwrites vector b.

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,			                      /* row index of the matrix utm */
		j;			                      /* col index of the matrix utm */    
 
	for (i=0 ; i < dim ; i++ )
	{  
		for(j=0; j<i ; j++) 
			b[i]=b[i]-utm[(i)*dim+(j)]*b[j]; 
		b[i] = b[i]/utm[(i)*dim+(i)];  
	}
}

/***************************************************************************************/

void chsol(double *vx, double* vb, double *mpd, int dim)
/*
   Purpose:
      Solve the following matrix equating using cholesky decomposition

		 R^t R x = b             (1)

	  where R is a upper triangular matrix.

   Input :  
	  vb,	y vector in (1)
	  mpd,	positive definite matrix R^t R in (1)
	  dim,	dimension of mpd

   Output: 
	  vx,	solution of (1)
	  mpd,	mpd will contain the cholesky decomposition
			of the original matrix at exit.

   Precondition:
	  1. The input matrix mpd is a positive definite matrix.
	  2. memory for vx should be allocated in the calling program.

   Function calls other than C or NR utilities:
      dpdch()
	  subfd()
	  subbk()

   Comments:
       - dtdsm is specialized to handle tri-diagonal matrix 
	     matrix multiplication more efficiently.
  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{  
	/* check if vx, vb, mpd are null*/
	if ( vx==nullptr|| vb==nullptr|| mpd==nullptr )
	{
		std::printf("Error, chsol: Input Error\nnull pointer \n");
		std::exit(-1);
	}
	/* cholesky decomposes mpd		*/
	std::memcpy(vx,vb,dim*sizeof(double));
	dpdch(mpd,dim);									/* mc = R^t R		*/
	subfd(vx,mpd,dim);								/* R^t R vx = Q^ vy	*/ 
	subbk(vx,mpd,dim);								/* vx = solution	*/  
}

/***************************************************************************************/

void dtdmv(double *v2, double * mtd, double * v1, int rdim, int cdim)
/*
   Purpose:
      dtdmm performes the following matrix multiplication (mm)
	  operation
            
			 v2 =  mtd * v1

	 where mtd,		a tridiagonal matrix, and
		   v1,v2,   vectors

   Input :
	  v1,    vector to be multiplied
	  mtd,   tridiagonal matrix
	  rdim,  row dimension of mtd
	  cdim,  col dimension of mtd

   Output:
      v2,    multiplication of matrix mtd and vector v1

   Function calls other than C or NR utilities:
      None.

   Comments:
       - dtdmm is specialized to handle tri-diagonal matrix 
	     matrix multiplication more efficiently.

  Author: Jaehoon Seol
  Date of last revision: 2/15/08
*/
{
	int i,													/* row index */
		j,												 /* column index */  
		left,								/* left end  of the i-th row */
		right;								/* right end of the i-th row */  

	if (v1 == nullptr || v2 == nullptr || mtd == nullptr)
	{
		std::printf("Error, dtdmm: Input Error \n\tnull pointer\n");
		std::exit(-1);
	}  
	/* calculate v2 = mtd * v1		*/
	std::memset(v2,0,rdim*sizeof(double)); 
	for (i=0 ; i < rdim ; i++ )
	{
		/* left and right help to skip operations over zero 
		    values in mtd */
		left 	= (i-2>=0? i-2: 0);      
		right 	= (i+2<=cdim-1? i+2: cdim-1); 
		for (j=left; j <= right; j++ ) 
			v2[i] += mtd[(i)*cdim+(j)]*v1[j];  
	} 
}

/***************************************************************************************/
double linearPoly(double x0, double y0, double x1, double y1, double xvalue)
/* Purpose:                                                            
      Linear polynomial function f passing through (x0,y0) and (x1,y1) 
                                                                       
       f(x0,y0,x1,y1,xvalue)                                           
                y1 - y0                                                 
       = y0 + -----------*(xvalue-x0)                                  
                x1 - x0    

   Input:                                                               
      x0,y0,    x & y coordinate of a point p0                         
      x1,y1,    x & y coordinate of a point p1                         
      xvalue,   point at which to evaluate f() 

   Output:
      f(x0,y0,x1,y1,xvalue)

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 7/4/08
*/
{
	double slope,                                   /* slope of the line */
	       rvalue ;                                      /* return value */ 

	/* evaluate the function         */
	slope = (y1-y0)/(x1-x0) ;
	/* check that slope != 0         */
 	if (std::fabs(slope)<std::pow(10.0,3.0)*DBL_EPSILON) 
		std::printf("linearPoly: Input Warning, slope ~ 0 \n"); 
	rvalue = y0 + slope*(xvalue-x0);
	return rvalue;
}

/***************************************************************************************/

double cubicPoly(double xleft, double xright, 
	   double ai, double bi, double ci, double di, double xvalue)
/* Purpose:                                                            
      Cubic polynomial function f defined by ai,bi,ci,di, i.e.,		   
       f(ai,bi,ci,di,xvalue)                                          
       =ai+bi*(xvalue-xleft)+ci*(xvalue-xleft)^2+di*(xvalue-xleft)^3   
      defined over [xleft, xright]      

   Input:                                                              
      xleft,       left end of the domain for f()                      
      xright,      right end of the domain for f()                     
      ai,bi,ci,di, coefficients of f()                                 
      xvalue,      point at which to evaluate f() 

   Output:
      f(ai,bi,ci,di,xvalue)

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 7/3/08
*/
{
	double tempX,                                    /* temp. variable   */
	       rvalue ;                                  /* return value     */
	/* confirm xleft<=xvalue<=xright */
    if ( xvalue < xleft || xvalue >xright)
	{
		std::printf("cubicPoly: Input Error, cubic poly. used outside domain\n");
		std::printf("xleft = %f, xright= %f xvalue = %f\n", xleft, xright, xvalue);
	}
	/* evaluate the function         */
	tempX  = xvalue-xleft;
	rvalue = ai+bi*tempX+ci*std::pow(tempX,2.0)+di*std::pow(tempX,3.0); 
	return rvalue;
}

/***************************************************************************************/

void setupT(bool edist, double * mt, int n, double * h)
/*
   Purpose: 
     This function sets up a positive definite, (n-1)x(n-1) tridiagonal matrix  
	 defined by

	    t     = 2*(h     + h    )/3,         (1)
         i,i        i-1     i
	    t      = t        = h   /3           (2)
		 i,i+1    i+1,i      i
	 Refer to page 179, Reinsch (1967).

   Input  :
     h =  [  h ,h , ... h     ]
              0  1       n-1       
     
	      where  h    = x    -  x   , i=0, 1,2,...,n-1
	              i      i+1     i
 
     n     = number of elements in the vector h.
	 edist = true,      if h is equi-distance vector.
	       = false,     otherwise

   Output :
      mt, The matrix defined by (1) and (2)

   Function calls other than C or NR utilities:
      None.

   References :
     1. C.H. Reinsch, Smoothing by spline functions, 1967, Num. math. 10, p177-183. 
   Comments:
     -The size of memory for mt should be n*sizeof(double) and
	  the memory should be allocated before calling setup_matrixT.

  Author: Jaehoon Seol
  Date of last revision: 2/18/08

*/
{
	int i,j,                            /* row     and column loop index */
		left,                           /* left  end of non-zero elememt */
		right;                          /* right end of non-zero element */
	double h0,                              /* equi-distance of h vector */
			x;                                          /* temp variable */ 


	/* check that mt is not null allocated             */
	if ( mt == nullptr )
	{
		std::printf("Error: setup_matrixT \n\tnull pointer\n");
		std::exit(-1);
	}

	/* initialize all elements in matrix mt(,) by 0.0  */
	std::memset(mt,0,(n-1)*(n-1)*sizeof(double));

	/* initialize all non-zero element in matrix mt(,) */ 
	h0=h[0];     /* will be used if h is equi-distance */
	for (i=0 ; i< n-1; i++ )
	{
		left  = (i-1>=0   ? i-1:0);
		right = (i+1<=n-2 ? i+1:n-2);
		for (j=left ; j<=right; j++)
		{
			if (edist == true )
			{
				if(i==j)			                 /* diagonal element */ 
					mt[(n-1)*(i) +(j)]=4.0*h0/3.0;   /* refer to (1)     */  
				else   
					mt[(n-1)*(i) +(j)]=h0/3.0; 
			}
			else
			{
				if(i==j)				             /* diagonal element */
				{
					x=h[i]+h[i+1] ;                  /* refer to (1)     */
					mt[(n-1)*(i) +(j)]=x*2.0/3.0;    /* refer to (1)     */ 
				} 
				else if ( j==i+1 )
				{
					mt[(n-1)*(i) +(j)]=h[i+1]/3.0;
				}
				else if ( j==i-1 )
				{
					mt[(n-1)*(i) +(j)]=h[i]/3.0;
				}
			}
		}
	}
}

/***************************************************************************************/

void setupQt(bool edist, double * mqt, int n, double * h)
/*
   Purpose: 
     This function sets up a positive definite, (n+1)x(n-1) tridiagonal matrix  
	 defined by

	    q     = -1/h     -1/h                 (1)
         i,i        i-1      i
	    q        = 1/h                        (2)
		  i+1,i       i
	    q        = 1/h                        (3)
		 i-1,i        i-1
	  Refer to page 179, Reinsch (1967).

   Input  :
     h =  [  h ,h , ... h     ]
              0  1       n-1       
     
	      where  h    = x    -  x   , i=0, 1,2,...,n-1
	              i      i+1     i
 
     n     = number of elements in the vector h.
	 edist = true,      if h is equi-distance vector.
	       = false,     otherwise

   Output :
      mt, The matrix defined by (1) and (2)

   Function calls other than C or NR utilities:
      None.

   References :
     1. C.H. Reinsch, Smoothing by spline functions, 1967, Num. math. 10, p177-183. 
   Comments:
     - The size of memory for mt should be n*sizeof(double) and
	   the memory should be allocated before calling setup_matrixT.

  Author: Jaehoon Seol
  Date of last revision: 2/18/08
*/
{
	int i,j,							    /* row and column loop index */ 
		right;							/* right end of non-zero element */ 
	double h0;							    /* equi-distance of h vector */


	/* check that mqt is not null allocated             */
	if ( mqt == nullptr )
	{
		std::printf("Error: setup_matrixT \n\tnull pointer");
		std::exit(-1);
	}

	/* initialize all elements in matrix mqt(,) by 0.0  */
	std::memset(mqt,0,(n+1)*(n-1)*sizeof(double));

	/* initialize all non-zero element in matrix mt(,) */ 
	h0=h[0];					   /* will be used if h is equi-distance */
	for (i=0 ; i< n-1; i++ )
	{
		right  = (i+2<n+1? i+2: n);			        /* right limit index */ 
		for (j=i ; j<=right; j++)
		{
			if (edist == true )
			{ 
				h0=h[0];
				if(i==j)							/* diagonal element  */ 
					mqt[(n+1)*(i) +(j)]=1.0/h0;          /* refer to (1) */  
				else if ( j==i+1 )			 /* element 1 away from diag */
					mqt[(n+1)*(i) +(j)]=-2.0/h0 ; 
				else if ( j==i+2 )			 /* element 2 away from diag */
					mqt[(n+1)*(i) +(j)]=1.0/h0; 
			}
			else
			{
				if(i==j)					         /* diagonal element */ 
					mqt[(n+1)*(i) +(j)]=1.0/h[i];        /* refer to (1) */  
				else if ( j==i+1 )			 /* element 1 away from diag */
					mqt[(n+1)*(i) +(j)]=-(1.0/h[i+1]+1.0/h[i]); 
				else if ( j==i+2 )			 /* element 2 away from diag */
					mqt[(n+1)*(i) +(j)]= 1.0/h[i+1]; 
			}
		}
	}
}

/***************************************************************************************/

int sspline(double *x, double *y, double *dyi, int num, double s, double *mc)
/*
   Purpose:  
     This function calculates the coefficients

	     ai, bi, ci, di, i=0,1,...n 

	 for the cubic spline 
	                                  2        3
	     f(x) = ai + bi(x-xi)+ci(x-xi)+di(x-xi),   x in [xi,xi+1)

	 The cubic spline f(x) minimizes
                         2
	          Int[ f''(x)  ]dx

	 among all functions f(x) satisfying
	                                   
                        f(xi)-yi      2              2
	          Sum( ----------------- )  <= s,  f in C [x0, xn] ---(1)
			               dyi

   Description of the algorithm:     
	 * start with p=0.0                                   * * T       T  2           * * (1) cholesky decomposition R  R of Q  D  Q + pT    * * (2) compute u from R^t Ru = Q^t y and 
	 * v = D Q u, e=v^t v                     * * (3) if e is greater than S                         *
	 * compute f = u^t T u and g = w^t w             *
	 * where R^T w = T u                     *
	 * replace p by p+(e-(Se)^(.5))/(f-p*g)          *
	 * restart with step (1).                        *
	 * Otherwise                     
	 * (4) compute a = y-D v, c= pu                       *
	 * compute b and d according to (8) and (9)       * Input  : 
     x,   describes x-coordinates of input data.
	      (n+1)x1 matrix, i.e.,
		     x = {x0,x1,x2,...,xn}'
     y,   describes y-coordinates of input data.
	      (n+1)x1 matrix, i.e.,
		     y = {y0,y1,y2,...,yn}'
     dyi, estimates the standard deviation of the ordinate yi
	      (n+1)x1 matrix
     num, number of input data 
	 s,   fidelity constant which controls the closeness of f(xi) to yi in (1)
          This constant is scaled by (x[n]-x[0]+1) within the function to 
		  give "s" more consistent meaning across different applications.

   Precondition:  
     mc is a one-dimensional array of size 4*n*sizeof(double).
	 Enough memory for mc should be allocated before calling 
	 ssspline

   Output :
     mc,   coefficient matrix of the smoothing cubic spline returns 
	 number of iteration used to compute p. If the return value is equal
	 to -1, that means the max. iteration number(25) has been used.

   Note :
	 Additional break statement "if (e<=s) break;" at the end of the loop 
	 has been added to keep the matrix
	                  T  2
	                 Q  D  Q + pT  
	 positive definite. Withough this condition, the matrix may not be 
	 positive definite in some cases in which case Cholesky decomposition
	 can not be used. If the user is not comfortable with this, he/she should
	 consider using QR factorization rather than Cholesky decomposition.

   References :
     1. C.H. Reinsch, Smoothing by spline functions, 1967, Num. Math. 10, p177-183. 
	 2. C.H. Reinsch, Smoothing by spline functions. II, 1970, Num. Math. 16, p451-454 

  Author: Jaehoon Seol
  Date of last revision: 8/10/08
*/
{
	 int	n,                                       /* col. dim of cmat */
			i,											   /* loop index */
			j;											   /* loop index */  
	 int	ii=0;									  /* main loop index */
	 double	e=0.0,                                            /* e=v^t v */  
			np=0.0,                                    /* new, updated p */
			p=1.0,                         /* p<-p+(e-(Se)^(.5))/(f-p*g) */
			g=0.0,										    /* g = w^t w */ 
			f=0.0;                                        /* f = u^t T u */
	 double *vh,                                    /* h_i = x_(i+1)-x_i */
			*va,	                            /* cubic spline coeff. a */
			*vc,	                            /* cubic spline coeff. c */
			*vu,	                               /* u: R^t R u = Q^t y */ 
			*vw,	                                 /* w: R^t R w = T u */ 
			*vtw,	                                          /* tw: T w */ 
			*vy2,	                                        /* y2: Q^t y */ 
			*vv,	                                         /* v= D Q u */
			*vtu,	                            /* temporary of size n-1 */ 
			*mq,	                                         /* matrix Q */
			*mqt,	                                       /* matrix Q^t */
			*mpt,	                                        /* mpt = p T */
			*mt,	                                         /* matrix T */
			*mpt1,	                         /* part 1:(fixed) Q^t D^2 Q */
			*mpt2,	                          /* part 2: Q^t D^2 Q + p T */ 
			*mtmp;	                                /* temp matrix for Q */ 
	 if (x==nullptr || y==nullptr || dyi==nullptr ||mc==nullptr)
	 {
		std::printf("Error, sspline: Input Error \n\tnull pointer\n");
		return 0;
	 }
	 n=num-1;	                          
	 s = s*(x[n]-x[0]+1);  

	 va		= static_cast<double *>(std::malloc((n+1)*sizeof(double)));
	 vc		= static_cast<double *>(std::malloc((n+1)*sizeof(double)));
	 vu		= static_cast<double *>(std::malloc((n-1)*sizeof(double)));
	 vw		= static_cast<double *>(std::malloc((n-1)*sizeof(double)));
	 vtw	= static_cast<double *>(std::malloc((n-1)*sizeof(double))); 
	 vh		= static_cast<double *>(std::malloc(n*sizeof(double)));
	 vv		= static_cast<double *>(std::malloc((n+1)*sizeof(double)));
	 vy2	= static_cast<double *>(std::malloc((n+1)*sizeof(double)));	    
	 vtu	= static_cast<double *>(std::malloc((n-1)*sizeof(double))); 	    
	 mq		= static_cast<double *>(std::malloc((n-1)*(n+1)*sizeof(double)));          
	 mt		= static_cast<double *>(std::malloc((n-1)*(n-1)*sizeof(double)));         
	 mpt	= static_cast<double *>(std::malloc((n-1)*(n-1)*sizeof(double)));
	 mqt	= static_cast<double *>(std::malloc((n-1)*(n+1)*sizeof(double)));        
	 mpt1	= static_cast<double *>(std::malloc((n-1)*(n-1)*sizeof(double)));  
	 mpt2	= static_cast<double *>(std::malloc((n-1)*(n-1)*sizeof(double)));  
	 mtmp	= static_cast<double *>(std::malloc((n-1)*(n+1)*sizeof(double)));  

	 /* setup vectors and matrices that does not *
	  * change thru. loop						 */
	 for(i=0; i<n; i++) vh[i] = x[i+1]-x[i];  
	 setupQt(false,mqt,n,vh);			            /* setup matrix Q^t  */ 
	 dtdmt(mq,mqt,n-1,n+1);				               /* setup matrix Q */
	 setupT(false,mt,n,vh);				               /* setup matrix T */ 
	 std::memset(mtmp,0,(n-1)*(n+1)*sizeof(double)); 
	 for(i=0; i<n+1; i++)				               /* compute  D^2 Q */ 
		 for(j=0; j<n-1; j++)
			 mtmp[(i)*(n-1)+(j)]=dyi[i]*dyi[i]*mq[(i)*(n-1)+(j)]; 
	 dtdmm(mpt1,mqt,mtmp,n-1,n+1,n-1);	             /* mpt1 = Q^t D^2 Q */  
	 std::free(mtmp); 
	 dtdmv(vy2,mqt,y,n-1,n+1);			           /* compute y2 = Q^t y */ 
 	 ii=0;
	 while(std::fabs(np-p)> std::pow(10.0,3.0)*DBL_EPSILON && ii++<35)
	 {      
		/* cholesky decomposition R^t R	of	*
		 * (mpt2=) Q^t D^2 Q + p*T		*/     
		p=np;
		dtdsm(mpt,p,mt,n-1,n-1);   
		dtdma(mpt2,mpt1,mpt,n-1,n-1);      
		chsol(vu,vy2,mpt2,n-1);		                 /* R^t R u = Q^t y	 */  
		dtdmv(vv,mq,vu,n+1,n-1);	                 /* compute v=D Q u	 */
		for(i=0; i<n+1 ; i++)	vv[i] = dyi[i]*vv[i]; 
		e=0.0 ; 
		for(i=0; i<n+1; i++)	e += vv[i]*vv[i];  
		dtdmv(vtu,mt,vu,n-1,n-1); 
		/* compute w in R^t R w = T u */
		std::memcpy(vw,vtu,(n-1)*sizeof(double));
		subfd(vw,mpt2,n-1);                                 /* vw has w  */ 
		subbk(vw,mpt2,n-1);                                 /* vw has w  */ 
		dtdmv(vtw,mt,vw,n-1,n-1); 
		/* compute f, g, and p		  */   
		f = 0.0; for(i=0; i<n-1; i++)	f += vu[i]*vtu[i];
		g = 0.0; for(i=0; i<n-1; i++)	g += vu[i]*vtw[i];
		if (s != 0.0 )
			np = p + std::sqrt(e/s)*(e-std::sqrt(s*e))/(f-p*g);
		else
			np = p + (e-std::sqrt(s*e))/(f-p*g);   
		if (e<=s) break;
	 } 
	 /* return matrix mc:						* * mc(0, ) : coefficients a				*
	  * mc(1, ) : coefficients b				*
	  * mc(2, ) : coefficients c				* * mc(3, ) : coefficients d				*/  
	 va[0]=y[0]-dyi[0]*vv[0]; 
	 vc[0]=0.0;	 
	 /* return coefficients b and d				*/ 
	 for (i=0; i < n ; i++)
	 {
		 va[i+1] = y[i+1]-dyi[i+1]*vv[i+1]; 
		 vc[i+1] = (i==n-1? 0:np*vu[i]);		
		 mc[(3)*(n)+(i)] = (vc[i+1]-vc[i])/(3.0*vh[i]) ;		  
		 mc[(1)*(n)+(i)] = (va[i+1]-va[i])/vh[i]-(vc[i]+mc[(3)*(n)+(i)]*vh[i])*vh[i];   
	 }  
	 /* return remaining coefficients			*/ 
	 std::memcpy(mc,va,n*sizeof(double));  
	 std::memcpy(mc+2*n,vc,n*sizeof(double)); 
	 /* free all working areas */
	 std::free(va);
	 std::free(vc);
	 std::free(vh);
	 std::free(vu);
	 std::free(vv);
	 std::free(vw);
	 std::free(vtw);
	 std::free(vtu);
	 std::free(vy2);
	 std::free(mq);
	 std::free(mt);
	 std::free(mqt); 
	 std::free(mpt);
	 std::free(mpt1); 
	 std::free(mpt2);  
	 if (ii==25) ii = -1; 
	 return ii;
}

/***************************************************************************************/

void postSmooth(double *xvalues, double *yvalues, double *dyi, int num1, double s, 
	 int xlow, int xhigh, double ky, double *vectX, int num2, double *vectY, double *cmat)
/* Purpose/functionality
...
*/
{
	int		i,j,							               /* loop index */
			numcoeff;					/* number of cubic spline pieces */
	double	dy_xlow,										/* dy(xlow)  */
			dy_xhigh,										/* dy(xhigh) */ 
			kx,          /* maximum x value (Kx in Kolen & Brennan, 2004)*/ 
			xvalue;					     /* temp. var. to store vectX[i] */  

	/* compute coefficient matrix of the cubic spline */
	numcoeff	= xhigh-xlow;  
    sspline(&xvalues[xlow],&yvalues[xlow],&dyi[xlow],numcoeff+1,s,cmat);
	dy_xlow		= cubicPoly(xvalues[xlow],xvalues[xlow+1],cmat[0],cmat[numcoeff],
			      cmat[2*numcoeff],cmat[3*numcoeff],xvalues[xlow]);
	dy_xhigh	= cubicPoly(xvalues[xhigh-1],xvalues[xhigh],cmat[numcoeff-1],cmat[2*numcoeff-1],
		          cmat[3*numcoeff-1],cmat[4*numcoeff-1],xvalues[xhigh]);
	kx			= xvalues[num1-1];		 
	for (i=0; i< num2; i++)
	{
		xvalue = vectX[i];
		/* Use Eq. 1,   */ 
		if ( -0.5 <= xvalue && xvalue < xvalues[xlow] ) 
			vectY[i] = linearPoly(-0.5,-0.5,xvalues[xlow],dy_xlow,xvalue);
		/* Use Eq. 3,   */ 
		else if ( xvalues[xhigh] < xvalue && xvalue <= kx+0.5 ) 
			vectY[i] = linearPoly(xvalues[xhigh],dy_xhigh,kx+0.5,ky+0.5,xvalue);
		/* Use Eq. 2,   */ 
		else if ( xvalues[xlow] <= xvalue && xvalue <= xvalues[xhigh] )
		{
			/* find the index of the  sub-interval     *
			 * to which xvalue belongs to              */
			for( j=0; j<numcoeff ;j++)
				if ( xvalues[xlow+j] <= xvalue && xvalue <= xvalues[xlow+j+1] ) 
					break;
			/* evaluate dy(xvalue)                     */ 
			vectY[i] = cubicPoly(xvalues[xlow+j],xvalues[xlow+j+1],cmat[j],cmat[j+numcoeff],
				       cmat[j+2*numcoeff],cmat[j+3*numcoeff],xvalue);
		}
		else
		{
			std::printf("postSmooth: Input Error\n");
			std::printf("xvalue = %8.5f xvalues[xlow]=%8.5f xvalues[xhigh]=%8.5f \n", 
				xvalue,xvalues[xlow],xvalues[xhigh]);
			std::exit(-1);
		}
	}
}

/***************************************************************************************/

double inverseCubicPoly(double left, double right, 
	   double ai, double bi, double ci, double di, double yvalue)
/* Purpose: ... */ 
{                  
	double diff,									   /* |right - left| */
	       side1,									   /* temp. variable */
		   side2,								       /* temp. variable */
		   left0=left,			   /* left end of cubic spline's domain  */
		   right0=right,           /* right end of cubic spline's domain */
		   mid = 0.0;                    /* middle point of left &right        */ 
	const double ERROR=std::pow(10.0,6.0)*DBL_EPSILON;
 
	/* check the precondition */
	side1 = cubicPoly(left,right,ai,bi,ci,di,left);
	side2 = cubicPoly(left,right,ai,bi,ci,di,right);
	if ( side1 >= side2 )
	{
		std::printf("inverseCubicPoly: Input Error 1\n");
		std::printf("left end = %f right end = %f\n", side1, side2);
		return yvalue;
	} 
	if ( yvalue < side1 || side2 <yvalue ) 
	{
		std::printf("inverseCubicPoly: Input Error 2\n");
		std::printf("left end=%f right end=%f yvalue=%f\n", side1, side2, yvalue);
		return yvalue;
	}  
	diff = std::fabs(right-left);    /* difference of the interval */ 
	while( diff > ERROR )
	{
		mid = (left+right)/2.0;
		side1 = cubicPoly(left0,right0,ai,bi,ci,di, mid)-yvalue;
		side2 = cubicPoly(left0,right0,ai,bi,ci,di, right)-yvalue;
		if (side1*side2 <= 0)
			left = mid;
		else
			right = mid;
		diff = std::fabs(right-left); 
	}
	return mid;
} 

/***************************************************************************************/

void inverseSSpline(double *ynodes, double * cmat, int n2,double *vectX, int n1,double *vectY)
/*
   Purpose: ...
 */
{
	int		 i,j;		                                /* loop index    */
	double   *xnodes;                                   /* x node values */

	/* setup y[i]=cubicPoly(,,,ynodes[i])         */
	xnodes=static_cast<double *>(std::malloc((n2+1)*sizeof(double)));
	for(j=0; j<n2; j++)
		xnodes[j]=cubicPoly(ynodes[j],ynodes[j+1],
			 cmat[j],cmat[n2+j],cmat[2*n2+j],cmat[3*n2+j],ynodes[j]);
	j--;
	xnodes[n2]=cubicPoly(ynodes[j],ynodes[j+1],
		       cmat[j],cmat[n2+j],cmat[2*n2+j],cmat[3*n2+j],ynodes[n2]);

	for (i=0; i< n1; i++)
	{
		/* evaluate the inverse of the cubic spline at vectY[i] */
		/* step 0: check xnodes[0]<= vectX[i]<= xnodes[n2]      */
		if ( vectX[i] <xnodes[0] || vectX[i] > xnodes[n2] )
		{
			std::printf("inverseSSpline: Input Error\n");
			std::printf("value = %f\n", vectX[i]);
			return;
		}

		/* step 1: find node index to which vectX[i] belongs    */ 
        for(j=0; j<n2; j++)  
			if (xnodes[j]<=vectX[i] && vectX[i]<=xnodes[j+1])
				break;	

		/* step 2: construct cubic spline using the coefficients *
		 * stored in cmat                                */  
		vectY[i]=inverseCubicPoly(ynodes[j],ynodes[j+1],
			      cmat[j],cmat[n2+j],cmat[2*n2+j],cmat[3*n2+j],vectX[i]);
	}
	std::free(xnodes);
}

/***************************************************************************************/

void inversePostSmooth(double *yvalues, double *xvalues, double *dxi, int num1, double s, 
	 int ylow, int yhigh, double kx, double *vectX, int num2, double *vectY, double *cmat)
/* Purpose/functionality ... */
{
	int		i,				                               /* loop index */
			numcoeff;                   /* number of cubic spline pieces */
	double	dx_ylow,		                                /* dx(ylow)  */
			dx_yhigh,		                                /* dx(yhigh) */ 
			ky,			 /* maximum y value (Ky in Kolen & Brennan, 2004)*/ 
			temp;                        /* temp. var. to store vectX[i] */  

	/* compute coefficient matrix of the cubic spline */

	numcoeff	= yhigh-ylow; 
    sspline(&yvalues[ylow],&xvalues[ylow],&dxi[ylow],numcoeff+1,s,cmat);
	dx_ylow		= cubicPoly(yvalues[ylow],yvalues[ylow+1],cmat[0],cmat[numcoeff],
			      cmat[2*numcoeff],cmat[3*numcoeff],yvalues[ylow]);
	dx_yhigh	= cubicPoly(yvalues[yhigh-1],yvalues[yhigh],cmat[numcoeff-1],cmat[2*numcoeff-1],
		          cmat[3*numcoeff-1],cmat[4*numcoeff-1],yvalues[yhigh]);
	ky			= yvalues[num1-1];

	for (i=0; i< num2; i++)
	{
		temp = vectX[i]; 
		if ( -0.5 <= temp && temp < dx_ylow )			/* inverse of Eq. 1 */
			vectY[i] = linearPoly(-0.5,-0.5,dx_ylow,yvalues[ylow],temp); 
		else if ( dx_yhigh < temp && temp <= kx+0.5 )   /* inverse of Eq. 3 */
			vectY[i] = linearPoly(dx_yhigh,yvalues[yhigh],kx+0.5,ky+0.5,temp); 
		else if ( dx_ylow <= temp && temp <= dx_yhigh ) /* inverse of Eq. 2 */		
			inverseSSpline(&yvalues[ylow],cmat,numcoeff,&vectX[i],1,&vectY[i]); 
		else
		{
			/* temp=vectX[i] should not be here */ 
			std::printf("inversePostSmooth: Input Error\n");
			std::printf("xvalue = %8.5f %8.5f %8.5f \n", temp,dx_ylow,dx_yhigh);
			std::exit(-1);
		} 
	}
}

/***************************************************************************************/