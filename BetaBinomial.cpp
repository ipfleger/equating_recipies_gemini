/*
  BetaBinomial.cpp
  Contains functions used for:
      (a) beta binomial smoothing--there are three types:
            2 parameter beta, binomial errors
            4 parameter beta, binomial errors
            4 parameter bete, compound binomial errorsin computing
      (b) equating beta-binomial smoothed x to scale of
          beta-binomial smoothed y

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

#include "ERutilities.h"
#include "NRutilities.h"
#include "BetaBinomial.h"
#include <cstdio>
#include <cstring>
#include <cmath>

/*************************************************************************/

void Wrapper_Smooth_BB(struct USTATS *x, int nparm, double rel,
	struct BB_SMOOTH *s)
	/*
	  Wrapper to do beta binomial smoothing.  There are three types:
		2 parameter beta, binomial errors
		4 parameter beta, binomial errors
		4 parameter bete, compound binomial errors

	  Input
		x     =  UTATS structure
		nparm = number of parameters for beta (2 or 4)
		rel   = reliability (usually KR20);
				if rel == 0, binomial errors;
				otherwise compund binomial errors

	  Output
		populates s, which is a BB_SMOOTH structure

	  Function calls other than C or NR utilities:
		Smooth_BB

	  R. L. Brennan

	  Date of last revision: 6/30/08

	*/
{
	Smooth_BB(x->n, x->ns - 1, x->dbl_fd, x->mts, nparm, rel, s);
}

/**************************************************************************/

void Smooth_BB(int n, int nitems, double *fd, double *mts,
	int nparm, double rel, struct BB_SMOOTH *s)
	/*
	  Performs beta binomial smoothing.  There are three types:
		2 parameter beta, binomial errors
		4 parameter beta, binomial errors
		4 parameter bete, compound binomial errors

	  With some exceptions, code follows procedures and equations in:
	  Hanson, B. A. (1991) Method of moments estimates for the four
	  parameter beta compund binomial model and the calculation of
	  classification consistency estimates. (ACT Research Report 91-5).
	  Iowa City, IA: American College Testing

	  Input
		n = number of persons
		items = number of items
		fd[] = frequency distribution
		mts[] = raw score moments
		nparm = number of parameters for beta (2 or 4)
		rel = reliability (usually KR20);
			  if rel == 0, binomial errors;
			  otherwise compund binomial errors

	  Output
		populates struct BB_SMOOTH s

	  Function calls other than C or NR utilities:
		Beta2Smooth()
		Beta4Smooth()
		CalcLordk()
		cum_rel_freqs()
		perc_rank()

	  R. L. Brennan and B. A. Hanson

	  Date of last revision: 6/30/08

	*/
{
	int i;

	s->num_persons = n;                                        /* # persons */
	s->num_items = nitems;                                       /* # items */
	for (i = 0; i <= 3; i++) s->rmoments[i] = mts[i];               /* raw moments */
	s->nparm = nparm;                              /* # parameters (2 or 4) */
	s->rel = rel;                     /* reliability --- almost always kr20 */
	s->density = dvector(0, nitems);                       /* fitted density */
	s->crfd = dvector(0, nitems);                /* fitted cum rel freq dist */
	s->prd = dvector(0, nitems);            /* fitted perc rank distribution */

	if (nparm == 2) {
		Beta2Smooth(fd, s);
	}
	else {
		s->lordk = (rel > 0) ? CalcLordk(rel, nitems, mts) : 0.;
		if (s->lordk < 0) s->lordk = 0;
		Beta4Smooth(fd, s);
	}

	cum_rel_freqs(0, static_cast<double>(nitems), 1, s->density, s->crfd);
	for (i = 0; i <= nitems; i++)
		s->prd[i] = perc_rank(0, static_cast<double>(nitems), 1, s->crfd, static_cast<double>(i));
}

/***********************************************************************/

void Print_BB(FILE *fp, char tt[], struct USTATS *x, struct BB_SMOOTH *s)
/*
  print results in BB_SMOOTH s
  
  Input
    fp   = file pointer for output
    tt[] = user supplied text identifier
    x    = struct USTATS 
    s    = struct BB_SMOOTH

  Function calls other than C or NR utilities: None

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
	int i;

	std::fprintf(fp, "\n\n%s\n\n", tt);
	std::fprintf(fp, "Input filename:  %s\n\n", x->fname);

	std::fprintf(fp, "Beta-Binomial Smoothing:  ");
	std::fprintf(fp, "%1d-parameter beta with ", s->nparm);
	std::fprintf(fp, "%s binomial error model", (s->lordk > 0) ? "compound" : "");
	std::fprintf(fp, "\n\n       Number of examinees: %6d", s->num_persons);
	std::fprintf(fp, "\n           Number of items: %6d", s->num_items);
	std::fprintf(fp, "\nReliability (usually KR20): %12.5f", s->rel);
	std::fprintf(fp, "\n                  Lord's k: %12.5f", s->lordk);
	std::fprintf(fp, "\n     Number of momemts fit: %6d", s->momentsfit);

	std::fprintf(fp, "\n\n***Parameter Estimates for Beta Distribution***");
	std::fprintf(fp, "\n\n      alpha: %12.5f", s->beta[0]);
	std::fprintf(fp, "\n       beta: %12.5f", s->beta[1]);
	std::fprintf(fp, "\nlower limit: %12.5f", s->beta[2]);
	std::fprintf(fp, "\nupper limit: %12.5f", s->beta[3]);

	std::fprintf(fp, "\n\n***Moments***");
	std::fprintf(fp, "\n\n              Raw      Fitted-Raw              True\n");
	std::fprintf(fp, "\nMean %12.5f    %12.5f      %12.5f", s->rmoments[0],
		s->fmoments[0], s->tmoments[0]);
	std::fprintf(fp, "\nS.D. %12.5f    %12.5f      %12.5f", s->rmoments[1],
		s->fmoments[1], s->tmoments[1]);
	std::fprintf(fp, "\nSkew %12.5f    %12.5f      %12.5f", s->rmoments[2],
		s->fmoments[2], s->tmoments[2]);
	std::fprintf(fp, "\nKurt %12.5f    %12.5f      %12.5f", s->rmoments[3],
		s->fmoments[3], s->tmoments[3]);

	std::fprintf(fp, "\n\n***Chi-Square Statistics for Fitted Distribution***");
	std::fprintf(fp, "\n\nLikelihood ratio: %12.5f", s->lrchisq);
	std::fprintf(fp, "\n         Pearson: %12.5f", s->pchisq);

	std::fprintf(fp, "\n\n***Frequencies and Proportions***");
	std::fprintf(fp, "\n\nScore   freq  fitted-freq          prop  fitted-prop"
		"  fitted-crfd   fitted-prd\n");
	for (i = 0; i <= s->num_items; i++) {
		std::fprintf(fp, "\n%5d %6.0f %12.5f       %7.5f %12.5f %12.5f %12.5f", i,
			x->dbl_fd[i], s->num_persons*s->density[i],
			x->rfd[i], s->density[i], s->crfd[i], s->prd[i]);
	}

	std::fprintf(fp, "\n\n");
	for (i = 1; i <= 78; i++) std::fprintf(fp, "*");
}

/************************************************************************/

void Wrapper_RB(char design, char method, char smoothing,
	struct USTATS *x, struct USTATS *y,
	struct BB_SMOOTH *bbx, struct BB_SMOOTH *bby, int rep,
	struct PDATA *inall, struct ERAW_RESULTS *r)
	/*
	  Wrapper for doing equipercentile equating with RG design
		and beta-binomial smoothing
		... (comments omitted for brevity) ...
	*/
{
	/* method name --- 10 characters; right justified */
	char *names[] = { "   y-equiv" };

	inall->rep = rep;               /* should be set to 0 for actual equating */
						/* counting of replications done in Wrapper_Bootstrap() */

	/* allocation and assignments for struct PDATA inall */

	if (inall->rep == 0) {     /* no assignment or stor alloc for bootstrap reps */
		std::strcpy(inall->xfname, x->fname);
		std::strcpy(inall->yfname, y->fname);
		inall->x = x;
		inall->y = y;
		inall->design = design;
		inall->method = method;
		inall->smoothing = smoothing;

		inall->nm = 1;
		inall->names = cmatrix(0, inall->nm - 1, 0, 11);    /* only one row/method, 0 */
		std::strcpy(inall->names[0], names[0]);

		inall->min = x->min;
		inall->max = x->max;
		inall->inc = x->inc;
		inall->fdx = x->fd;
		inall->n = x->n;

		inall->bbx = bbx;
		inall->bby = bby;
	}

	/* allocation and assignments for r */

	if (inall->rep <= 1) {         /* no storage allocation for bootstrap reps >1 */
		r->eraw = dmatrix(0, 0, 0, loc(inall->max, inall->min, inall->inc));
		r->mts = dmatrix(0, 0, 0, 3);                    /* 0,3 is for the 4 moments */
	}

	/* Compute equating results */

	EquiEquate(y->ns, y->min, y->inc, bby->crfd, x->ns, bbx->prd, r->eraw[0]);

	/* get moments */

	MomentsFromFD(inall->min, inall->max, inall->inc, r->eraw[0], inall->fdx, r->mts[0]);

}

/********************************************************************************/

void Print_RB(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
	/*
	  print RB results (RG; equipercentile; beta-binomial smoothing)
	  ... (comments omitted for brevity) ...
	*/
	int i, j;

	std::fprintf(fp, "\n\n%s\n\n", tt);

	if (inall->rep > 0)  std::fprintf(fp, "Bootstrap relication %d\n\n", inall->rep);

	std::fprintf(fp, "Equipercentile Equating with Random Groups Design:");

	std::fprintf(fp, "\n\nBeta-Binomial Smoothing for x (new form): ");
	std::fprintf(fp, "\n   %1d-parameter beta with ", inall->bbx->nparm);
	std::fprintf(fp, "%s binomial error model", (inall->bbx->lordk > 0) ? "compound" : "");
	std::fprintf(fp, "\nBeta-Binomial Smoothing for y (old form): ");
	std::fprintf(fp, "\n   %1d-parameter beta with ", inall->bby->nparm);
	std::fprintf(fp, "%s binomial error model", (inall->bby->lordk > 0) ? "compound" : "");

	std::fprintf(fp, "\n\nInput file for x: %s\n", inall->xfname);
	std::fprintf(fp, "Input file for y: %s\n\n", inall->yfname);

	for (i = 1; i <= 30; i++) std::fprintf(fp, "-");

	/* following code set up for any number of methods */

	std::fprintf(fp, "\n\n");
	for (j = 1; j <= 15 + (inall->nm * 12 - 18) / 2; j++) std::fprintf(fp, " ");
	std::fprintf(fp, "Equated Raw Scores");
	std::fprintf(fp, "\n\n               ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "   Method %d:", j);
	std::fprintf(fp, "\nRaw Score (x)  ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "  %s", inall->names[j]);
	std::fprintf(fp, "\n\n");

	for (i = 0; i <= loc(inall->max, inall->min, inall->inc); i++) {
		std::fprintf(fp, "\n %12.5f  ", score(i, inall->min, inall->inc));
		for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f", r->eraw[j][i]);
	}

	std::fprintf(fp, "\n\n         Mean  ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f", r->mts[j][0]);
	std::fprintf(fp, "\n         S.D.  ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f", r->mts[j][1]);
	std::fprintf(fp, "\n         Skew  ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f", r->mts[j][2]);
	std::fprintf(fp, "\n         Kurt  ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f", r->mts[j][3]);

	std::fprintf(fp, "\n\n");
	for (j = 1; j <= 63; j++) std::fprintf(fp, "*");
}

/***********************************************************************/

short Beta4Smooth(double *rfreq, struct BB_SMOOTH *betafit)
/*
	Calculate 4 parameter beta (compound) binomial distribution.

	Input
	
		rfreq			         = Raw (# correct) score frequencies
		betafit->num_items = Number of items
		betafit->lordk		 = Lord's k
	
	Output
		betafit->density	 = fitted raw score distribution (density)
                         (space allocated before function called)

  Returns zero if no error, nonzero if error.

  Function calls other than C or NR utilities:
    BetaMoments() 
    CalcBetaPara()
    ObsDensity()
    ObsDenK()
    LRChiSqr()
    PChiSqr()
    MomentsFromRFD()
                                                
   B. A. Hanson with updates by R. L. Brennan 

  Date of last revision: 6/30/08  
*/
{
	short tfit=-1;		       /* number of moments fit to obtain estimates */
    int nsamp;                                           /* number of persons */
	double *sdist=nullptr,		             /* pointer to smoothed proportions */
		   *sfreq=nullptr;		         /* pointer to smoothed frequencies */
    double nctmoments[4];	                  /* non-central true score moments */
	int nitems,		                             /* number of items on test */
		  s,			                                      /* loop index */
	    nlrchi2, npchi2;	   /* number of categories used for chi-squares */
	
    nsamp = betafit->num_persons;
	nitems = betafit->num_items;
	sdist = betafit->density;	
    sfreq = dvector(0,nitems);                                  /* allocation */	

    /* calculate true score moments */

	BetaMoments(nitems,betafit->lordk,betafit->rmoments,betafit->tmoments,
              nctmoments);

    /* calculate parameters of beta true score distribution  */

	tfit = CalcBetaPara(nitems,betafit->tmoments,nctmoments,betafit->beta);
	
    /* calculate observed score density */

	if(ObsDensity(nitems+1,nsamp,betafit->beta,sfreq)) 
	{
		/* Assign uniform distribution and exit */
		for (s=0; s<=nitems; s++) sdist[s] = 1.0 / static_cast<double>(nitems+1);
		tfit = 0;
		goto EXIT;
	}

    /* adjust for compound binomial if k != 0 */

	if (betafit->lordk != 0.0) ObsDenK(nitems,betafit->lordk,sfreq);

	/* compute smoothed proportions */

	for (s=0; s<=nitems; s++) sdist[s] = sfreq[s] / nsamp;

    /* compute chi-square values */
	
	betafit->lrchisq = LRChiSqr(nitems+1,rfreq,sfreq,&nlrchi2);
	betafit->pchisq = PChiSqr(nitems+1,0.0,rfreq,sfreq,&npchi2);

    /* calculate fitted observed score moments */

    MomentsFromRFD(0, static_cast<double>(nitems),1,nullptr,sdist,betafit->fmoments);

EXIT:

	if (sfreq) free_dvector(sfreq,0,nitems);           /* deallocate  */

	betafit->momentsfit = tfit;
	
	if (tfit < 0) return -1;
	else return 0;
}

/*********************************************************************/

short Beta2Smooth(double *rfreq, struct BB_SMOOTH *betafit)
/*
	Calculate 2 parameter beta binomial distribution.
    ... (comments omitted) ...
*/
{
	short tfit=-1;		        /* number of moments fit to obtain estimates */
	int nitems,		                              /* number of items on test */
		  s,			                                       /* loop index */
	    nlrchi2, npchi2,	    /* number of categories used for chi-squares */
      nsamp;                                            /* number of persons */
	double *sdist=nullptr,		              /* pointer to smoothed proportions */
		   *sfreq=nullptr,		          /* pointer to smoothed frequencies */
	       *para;				     /* pointer to model parameter estimates */
	double nctmoments[4];	               /* non-central true score moments */
	

    nsamp = betafit->num_persons;
	nitems = betafit->num_items;
	sdist = betafit->density;
    sfreq = dvector(0,nitems);                                   /* allocation */
	betafit->lordk = 0.0;	                 /* binomial error distribution */		

    /* calculate true score moments */

	BetaMoments(nitems,betafit->lordk,betafit->rmoments,betafit->tmoments,
              nctmoments);

    /* calculate parameters of beta true score distribution  */

	para = betafit->beta;
	para[2] = 0.0;  para[3] = 1.0;
	if (!EstNegHypGeo(nitems,betafit->tmoments,para))
	{
		/* if two moments not fit, set alpha=1.0 and fit mean */
		para[2] = 0.0;  para[3] = 1.0;
		para[0] = 1.0;
		para[1] = (static_cast<double>(nitems) / (betafit->rmoments)[0]) - para[0];
		if (para[1] > 0.0)
		{
			tfit = 1;
		}
		else
		{
			/* if mean not fit return uniform distribution on 0 to 1.0 */
			para[1] = 1.0;
			tfit = 0;
		}
		
	}
	else tfit = 2;
	
    /* calculate observed score density */

	if(ObsDensity(nitems+1,nsamp,betafit->beta,sfreq)) 
	{
		/* Assign uniform distribution and exit */
		for (s=0; s<=nitems; s++) sdist[s] = 1.0 / static_cast<double>(nitems+1);
		tfit = 0;
		goto EXIT;
	}

	/* compute smoothed proportions */

	for (s=0; s<=nitems; s++) sdist[s] = sfreq[s] / nsamp;

    /* compute chi-square values */

	betafit->lrchisq = LRChiSqr(nitems+1,rfreq,sfreq,&nlrchi2);
	betafit->pchisq = PChiSqr(nitems+1,0.0,rfreq,sfreq,&npchi2);

    /* calculate fitted observed score moments */

    MomentsFromRFD(0,static_cast<double>(nitems),1,nullptr,sdist,betafit->fmoments);


EXIT:

	if (sfreq) free_dvector(sfreq,0,nitems);                 /* deallocate  */

	betafit->momentsfit = tfit;
	
	if (tfit < 0) return -1;
	else return 0;
}

/*****************************************************************************/

short CalcBetaPara(int nitems, double *moment, double *nctmoment, double *para)
/*
  Estimate parameters of 4-parameter beta binomial model
  ... (comments omitted) ...
*/
{
	
	/* Only try to fit more than one moment if the true score
		  s.d. is greater than zero */
	if (moment[1] > 0.0)
	{
	  /* try to fit all four moments */
		if (CalcBetaParaMM(nitems,moment,para)) return(4);
	
	  /* if four moments not fit try to fit three moments 
	     such that the squared difference in the estimated 
	     and observed kurtosis is as small as possible */
		if (CalcBetaParaLS(nitems,moment,nctmoment,para)) return(3);
	
	  /* if three moments not fit set lower = 0.0 and upper = 1.0 
	     and fit first two moments */
	 	para[2] = 0.0;  para[3] = 1.0;
		if (EstNegHypGeo(nitems,moment,para)) return(2);
	}
	
    /* if two moments not fit, set alpha=1.0 and fit mean */
	para[2] = 0.0;  para[3] = 1.0;
	para[0] = 1.0;
	para[1] = (static_cast<double>(nitems) / moment[0]) - para[0];
	if (para[1] > 0.0) return(1);

    /* if mean not fit return uniform distribution on 0 to 1.0 */
	para[1] = 1.0;
	return(0);
}

/***************************************************************************/

short CalcBetaParaMM(int nitems, double *moment, double *para)
/*
   Estimates parameters of 4-parameter beta binomial model
   by method of moments.
   ... (comments omitted) ...
*/
{
	double r,rad;
	double b1,b2,bma;	                  /* used to hold intermediate values */
	double l,u,a,b;		                        /* temporarily hold parameter */
	
    /* calculate alpha and beta */
	b2 = moment[3];
	b1 = moment[2] * moment[2];
	
	r = 6.0 * (b2-b1-1.0);
	r /= 6.0 + 3.0*b1 - 2.0*b2;                               /* calculate r */
	
	rad = (r+2.0) * (r+3.0) * b2;
	rad -= 3.0 * (r + 1.0) * (r - 6.0);
	rad = (24.0 * (r+1.0)) / rad;
	if (rad > 1.0) return(0);	             /* cannot compute alpha and beta */
	rad = std::sqrt(1.0-rad);
	
	r /= 2.0;
	if (moment[2] < 0.0)  	              /* compute values of alpha and beta */
	{
		a = r * (1.0 + rad);
		b = r * (1.0 - rad);
	}
	else
	{
		a = r * (1.0 - rad);
		b = r * (1.0 + rad);
	}
	
	if (a <= 0.0 || b <= 0.0) return(0);

    /* compute values of lower and upper limit parameters */
	bma = (a+b) * std::sqrt(a+b+1);
	bma *= moment[1];
	bma /= std::sqrt(a*b);
	
	l = -bma * (a / (a+b));
	l += moment[0];			                            /* compute lower limit */
	
	u = bma + l;			                            /* compute upper limit */
	
    /* assign values for output */
	para[2] = l / static_cast<double>(nitems);  
	para[3]= u / static_cast<double>(nitems); 
	para[0]=a; 
	para[1]=b;
	
	if (para[2] < 0.0 || para[3] > 1.0) return(0);
	
	return(1);
}

/********************************************************************************/

short CalcBetaParaLS(int nitems, double *tmoment, double *nctmoment, double *para)
/*
  Find value of lower limit of 4-parameter beta compound binomial model
  that minimizes squared difference between predicted and observed kurtosis,
  such that observed and predicted mean, variance and skewness are equal.
  ... (comments omitted for brevity) ...
*/
{
	double 	kurt,	/* squared difference of fitted and observed kurtosis */
			ll,tll = 0.0;		                            /* holds lower limits */
	double tpara[4];	           /* temporary space for beta parameters */
	double delta,rkurt,lkurt;
	int existr, existl, i;


		/* Initialize squared kurtosis difference to be infinity.
		   HUGE_VAL is a macro in <cmath> */

	kurt = HUGE_VAL;

	/* find solution for lower limit of 0 */

	existl = Kurtfuncd(0.0,nitems,para,tmoment,nctmoment,&lkurt);
	if (existl) kurt = lkurt;
	
	/* find solution for upper limit of 1 */

	existr = KurtfuncUpper(1.0,nitems,tpara,tmoment,nctmoment,&rkurt);
	if (existr && rkurt < kurt)
	{
		kurt = rkurt;
		for (i=0; i<4; i++) para[i] = tpara[i];
	}

	/* attempt to find a solution with a lower squared difference
		 in kurtosis than two solutions found above */

	delta = .05;
	ll = delta;
	existl = 0;		 /* existl will flag whether a better solution is found */
	while (ll < .55)
	{
		if (Kurtfuncd(ll,nitems,tpara,tmoment,nctmoment,&lkurt))
		{
			if (lkurt < kurt) /* a better solution is found if kurtl < kurt */
			{
				kurt = lkurt;
				for (i=0; i<4; i++) para[i] = tpara[i];
				tll = ll;
				existl = 1;
			}
		}				
		ll += delta;
	}
	
	
	if (kurt == HUGE_VAL) return(0);	/* no valid solution has been found */
	
	/* if no better solution than lower limit = 0 or
		 upper limit = 1 has been found then return */

	if (!existl) return(1);
	
	/* loop to find solution to accuracy of KACC */

	delta = .01;
	ll = tll;
	while (delta >= KACC)
	{
	  /* evaluate function at points to left and right 
       of current solution "ll" */

		existr = Kurtfuncd(ll+delta,nitems,para,tmoment,nctmoment,&rkurt);
		existl = Kurtfuncd(ll-delta,nitems,para,tmoment,nctmoment,&lkurt);
		if (existr && (rkurt < kurt))	               /* step to the right */
		{
			ll += delta;
			while(((ll+delta) < 1.0) &&
				Kurtfuncd(ll+delta,nitems,para,tmoment,nctmoment,&rkurt) && 
				(rkurt < kurt))
			{
				ll += delta;
				kurt = rkurt;
			}
		}
		else if (existl && (lkurt < kurt))	         /* step to the left */
		{
			ll -= delta;
			while(((ll-delta) > 0.0) &&
				Kurtfuncd(ll-delta,nitems,para,tmoment,nctmoment,&lkurt) && 
				(lkurt < kurt))
			{
				ll -= delta;
				kurt = lkurt;
			}
		}
		delta /= 10.0;
	}                                                /* end while loop */
				
	/* calculate final values to return */

	if (Kurtfuncd(ll,nitems,para,tmoment,nctmoment,&kurt)) return(1);

	return(0);	             /* this statement should never be reached */

}

/***********************************************************************/

short Kurtfuncd(double x, int nitems, double *para, double *tmoment,
	              double *nctmoment, double *kurt)
/*
  ... (comments omitted) ...
*/
{

	double tkurt;

	*kurt = -1.0;	              /* initialize in case of return of 0 */
	para[2] = x;
	if (!FindUpper(para,nctmoment)) return(0);
	if (!EstNegHypGeo(nitems,tmoment,para)) return(0);
	tkurt = CalcKurt(para);
	tkurt = tmoment[3] - tkurt;
	tkurt *= tkurt;
	*kurt = tkurt;
		
	return(1);	
}

/*********************************************************************/

short KurtfuncUpper(double x, int nitems, double *para, double *tmoment,
	double *nctmoment, double *kurt)
/*
   ... (comments omitted) ...
*/
{
	double tkurt;
	
	*kurt = -1.0;	              /* initialize in case of return of 0 */
	para[3] = x;
	if (!FindLower(para,nctmoment)) return(0);
	if (!EstNegHypGeo(nitems,tmoment,para)) return(0);
	tkurt = CalcKurt(para);
	tkurt = tmoment[3] - tkurt;
	tkurt *= tkurt;
	*kurt = tkurt;
		
	return(1);	
}

/**********************************************************************/

short FindUpper(double *para, double *tmoment)
/*
  ... (comments omitted) ...
*/
{
	double	fr1,	                             /* temporary */
						m1,m2,m3,  /* first through third central moments */
						lower;  /* lower limit of true score distribution */
	double	numerator;	       /* numerator of expression for upper limit */

	lower = para[2];
	m1 = tmoment[0];
	m2 = tmoment[1];
	m3 = tmoment[2];
	
	/* calculate upper limit */
	
	/* numerator */
	fr1 = m1 * m1 * (m2*lower - 2.0 * m3);
	fr1 += m2*m2 * (m1 - 2.0 * lower);
	fr1 += m3 * (m2 + m1 * lower);
	numerator = fr1;
	
	/* denominator */
	fr1 = m1*m1 * (2.0 * m1 * lower - m2);
	fr1 += m2 * (2.0*m2 - 3.0*m1*lower);
	fr1 += m3 * (lower - m1);
	
	para[3] = numerator / fr1;
	
	if (para[3] <= para[2] || para[3] > 1.0) return(0);	     /* invalid */
	
	return(1);
}

/************************************************************************/

short FindLower(double *para, double *tmoment)
/*
  ... (comments omitted) ...
*/
{
	double	fr1,	                             /* temporary */
						m1,m2,m3,  /* first through third central moments */
						upper;	/* upper limit of true score distribution */
	double	numerator;	       /* numerator of expression for upper limit */

	upper = para[3];
	m1 = tmoment[0];
	m2 = tmoment[1];
	m3 = tmoment[2];
	
  /* calculate lower limit */
	
	/* numerator */
	fr1 = m1 * m1 * (m2*upper - 2.0 * m3);
	fr1 += m2*m2 * (m1 - 2.0 * upper);
	fr1 += m3 * (m2 + m1 * upper);
	numerator = fr1;
	
	/* denominator */
	fr1 = m1*m1 * (2.0 * m1 * upper - m2);
	fr1 += m2 * (2.0*m2 - 3.0*m1*upper);
	fr1 += m3 * (upper - m1);
	
	para[2] = numerator / fr1;
	
	if (para[3] <= para[2] || para[2] < 0.0) return(0);  /* invalid */
	
	return(1);
}

/********************************************************************/

short EstNegHypGeo(int nitems, double *moment, double *para)
/*
  ... (comments omitted) ...
*/
{
	double smean, svar;	          /* standardized mean and variance */
	double maxmin;		  /* difference between minimum and maximum */
	double dn;
	double *ppara;
	
	dn = static_cast<double>(nitems);
	maxmin = dn * (para[3] - para[2]);
	smean = (moment[0] - dn*para[2]) / maxmin;
	svar = (moment[1]*moment[1]) / (maxmin*maxmin);
	
	ppara = para;
	*ppara = smean*smean * (1.0 - smean);
	*ppara /= svar;
	*ppara -= smean;
	
	ppara = para+1;
	*ppara = smean * (1.0-smean);
	*ppara /= svar;
	*ppara -= 1.0;
	*ppara -= para[0];
	if (para[0] > 0.0 && para[1] > 0.0)
		return(1);
	else
		return(0);
}
	
/*****************************************************************/

double CalcKurt(double *para)
/*
  Returns value of kurtosis give parameters alpha and beta 
  of 4 parameter beta distribution.
  ... (comments omitted) ...
*/
{
	double a,b;	                     /* beta parameters */
	double		                          /* components of kurtosis */
		k1, k2, k3, k4, k5, k6;
	double kurt;
		
    /* compute components used in calculations */
	a = para[0];
	b = para[1];
	k1 = 3.0 * (a+b+1.0);
	k2 = 2.0 * (a+b) * (a+b);
	k3 = a*b;
	k4 = a+b-6.0;
	k5 = a + b + 2.0;
	k6 = a+b+3.0;
		
    /* compute kurtosis */
	kurt = (k1 * (k2 + k3*k4)) / (k3*k5*k6);
	
	return(kurt);	
}

/*****************************************************************/

void BetaMoments(int n, double k, const double *rmoment, 
                 double *tmoment, double *nctmoment)
/*
  Given raw score mean, s.d., skewness and kurtosis,
  compute true score mean, s.d., skewness and kurtosis
  ... (comments omitted) ...
*/
{
	double 
		cmoment[4],		                          /* central moments */
		ncmoment[4],	                      /* non-central moments */
		fmoment[4];		                         /* factoral moments */
	double 				                   /* hold temporary results */
		prod2, prod3, n2, dn, dnm2, kr;

    /* compute raw score central moments */
	cmoment[2] = rmoment[2] * rmoment[1] * rmoment[1] * rmoment[1];
	cmoment[1] = rmoment[1] * rmoment[1];		         /* variance */
	cmoment[3] = rmoment[3] * cmoment[1] * cmoment[1];
	
    /* compute raw score non-central moments */
	prod2 = rmoment[0] * rmoment[0];
	ncmoment[1] = cmoment[1] + prod2;
	prod3 = prod2 * rmoment[0];
	ncmoment[2] = cmoment[2] + 3.0 * cmoment[1] * rmoment[0] + prod3;
	ncmoment[3] = cmoment[3] + 4.0 * rmoment[0] * cmoment[2] +
					6.0 * prod2 * cmoment[1] + prod3 * rmoment[0];

    /* compute raw factoral moments */
	fmoment[1] = ncmoment[1] - rmoment[0];
	fmoment[2] = ncmoment[2] - 3.0 * ncmoment[1] + 2.0 * rmoment[0];
	fmoment[3] = ncmoment[3] - 6.0*ncmoment[2] + 11.0*ncmoment[1] -
					6.0*rmoment[0];

    /* compute true score non-central moments */
	dn = static_cast<double>(n);
	n2 = static_cast<double>(n) * static_cast<double>(n-1);
	/* second moment */
	dnm2 = 1.0;
	kr = k * 2.0;
	ncmoment[0] = rmoment[0] / dn;
	nctmoment[0] = ncmoment[0];
	ncmoment[1] = (fmoment[1]/dnm2 + kr*ncmoment[0]) / (n2+kr);
	nctmoment[1] = ncmoment[1];
	/* third moment */
	dnm2 = dn-2.0;
	kr = k * 6.0;
	ncmoment[2] = (fmoment[2]/dnm2 + kr*ncmoment[1]) / (n2+kr);
	nctmoment[2] = ncmoment[2];
	/* fourth moment */
	dnm2 *= dn-3.0;
	kr = k * 12.0;
	ncmoment[3] = (fmoment[3]/dnm2 + kr*ncmoment[2]) / (n2+kr);
	nctmoment[3] = ncmoment[3];
	
	/* put true score moments on observed score scale */
	dn *= dn;
	ncmoment[1] *= dn;
	dn *= static_cast<double>(n);
	ncmoment[2] *= dn;
	dn *= static_cast<double>(n);
	ncmoment[3] *= dn;

    /* compute true score central moments */
	prod2 = rmoment[0] * rmoment[0];
	cmoment[1] = ncmoment[1] - prod2;
	prod3 = prod2 * rmoment[0];
	cmoment[2] = ncmoment[2] - 3.0 * ncmoment[1] * rmoment[0] + 2.0*prod3;
	cmoment[3] = ncmoment[3] - 4.0 * rmoment[0] * ncmoment[2] +
					6.0 * prod2 * ncmoment[1] - 3.0 * prod3 * rmoment[0];

    /* compute true score mean, s.d., skewness and kurtosis */
	tmoment[0] = rmoment[0];	
	/* The central true score moment may be negative, if so assign
	   the true score s.d., skewness and kurtosis to zero. */
	if (cmoment[1] > 0.0)
	{
		tmoment[1] = std::sqrt(cmoment[1]);
		tmoment[2] = cmoment[2] / (tmoment[1] * tmoment[1] * tmoment[1]);
		tmoment[3] = cmoment[3] / (cmoment[1] * cmoment[1]);
	}
	else
	{
		tmoment[1] = 0.0;
		tmoment[2] = 0.0;
		tmoment[3] = 0.0;
	}
}

/*****************************************************************/

short ObsDensity(int n, int nsamp, double *beta, double *scounts)
/*
  Calculate observed density given parameters of 4-parameter
  beta-binomial distribution. 
  ... (comments omitted) ...
*/
{
	double		              /* matrices used to fit distribution */
		*m1=nullptr,
		*m2=nullptr,
		**m3=nullptr,
		*m4=nullptr,
		*m5=nullptr,
		*prod1=nullptr,
		*prod2=nullptr;
	double *pm1, *pm2, *pm3;
	int i;	                                  /* loop index */
	double *ps;		                          /* pointer to scounts */
	double 	diff,	 /* difference in high and low cutoffs for beta */
			dconst;	    /* constant used to adjust computed density */
	int s;			                                  /* loop index */
	int m;		   /* temporary variable used in 2rd matrix product */
	short err=-1;		                          /* value returned */

	diff = beta[3] - beta[2];
	dconst = std::pow(diff, static_cast<double>(n - 1));
	if (dconst <= 0.0) return -2;	         /* check for underflow */
	dconst *= nsamp;

    /* allocate space for matrices */
	m1 = static_cast<double *>(std::malloc(n * sizeof(double)));
	if (!m1) goto EXIT;
	m2 = static_cast<double *>(std::malloc(n * sizeof(double)));
	if (!m2) goto EXIT;
	m4 = static_cast<double *>(std::malloc(n * sizeof(double)));
	if (!m4) goto EXIT;
	m5 = static_cast<double *>(std::malloc(n * sizeof(double)));
	if (!m5) goto EXIT;
	prod1 = static_cast<double *>(std::malloc(n * sizeof(double)));
	if (!prod1) goto EXIT;
	prod2 = static_cast<double *>(std::malloc(n * sizeof(double)));
	if (!prod2) goto EXIT;
	m3 = static_cast<double **>(std::malloc(n * sizeof(double *)));
	if (!m3) goto EXIT;
	for (s=0; s<n; s++) m3[s] = nullptr;	            /* initialize */
	for (s=0; s<n; s++)
	{
		m3[s] = static_cast<double *>(std::malloc((n-s) * sizeof(double)));
		if (!m3[s]) goto EXIT;
	}
	
	/* set err to zero to indicate memory successfully allocated */
	err = 0;
	
    /* calculate m1 */
	CalcM15(n,beta[0],m1);

    /* calculate m5 */
	CalcM15(n,beta[1],m5);

    /* calculate m2 */
	if ((err = CalcM24(n,beta[2]/diff,m2))) goto EXIT;

    /* calculate m4 */
	if ((err = CalcM24(n,(1.0-beta[3])/diff,m4))) goto EXIT;

    /* calculate M3 */
	CalcM3(n,beta,m3);

    /* loop to calculate beta-binomial frequencies */
	ps = scounts;
	for (s=0; s<n; s++)
	{
	    /* calculate first matrix product - M1 * M2 */
		i = s+1;
		pm1 = m1+n-1-s;
		pm2 = m2;
		pm3 = prod1;
		while (i--)
		{
			*pm3++ = *pm1++ * *pm2++;
		}

	    /* calculate second matrix product-- result of first product * M3 */
		pm1 = prod2;
		m = s + 1;
		for (i=0; i<n-s; i++)
		{
			*pm1++ = dotprod(m,prod1,m3[i]);
		}

  	    /* calculate third matrix product-- result of second product * M4 */
		pm1 = prod1;
		pm2 = prod2;
		pm3 = m4;
		i = n-s;
		while(i--)
		{
			*pm1++ = *pm2++ * *pm3++;
		}
		
	    /* calculate fourth matrix product-- result of third product * M5 */
		*ps++ = dconst * dotprod(n-s,prod1,m5+s);
		
	}                    /* end of loop for frequency calculation */	

/* release temporary space */
EXIT:
	if (m1) std::free(m1);
	if (m2) std::free(m2);
	if (m4) std::free(m4);
	if (m5) std::free(m5);
	if (prod1) std::free(prod1);
	if (prod2) std::free(prod2);
	if (m3)
	{
		for (s=0; s<n; s++)
		{
			if (m3[s]) std::free(m3[s]);
		}
		std::free(m3);
	}
	
	return err;
}

/*****************************************************************/

void ObsDenK(int n, double k, double *scounts)
/*
  ... (comments omitted) ...
*/
{
	double *pscounts;
	double d2,pxm1,px,pxp1,dsp1,kn2,dn;
	int s;
	
	dn = static_cast<double>(n);
	kn2 = k / (dn * static_cast<double>(n-1));
	
	/* raw score of 0 */
	pxp1 = (dn-1.0) * scounts[1];
	scounts[0] -= kn2 * pxp1;
	px = pxp1;
	pxm1 = 0.0;
	
	/* raw scores 1 through n-1 */
	pscounts = scounts+1;
	for (s=1; s<n; s++)
	{
		dsp1 = static_cast<double>(s+1);
		pxp1 = (dsp1) * (dn-dsp1) * *(pscounts+1);
		d2 = pxm1 - 2.0*px + pxp1;
		*pscounts -= kn2 * d2;
		pxm1 = px; px = pxp1; 
		pscounts++;
	}
	
	/* raw score n */
	scounts[n] -= kn2 * pxm1;

}

/**********************************************************************/

double dotprod(int n, double *v1, double *v2)
/*
  ... (comments omitted) ...
*/
{
	double prod;
	
	prod = 0.0;
	while(n--)
	{
		prod += *v1++ * *v2++;
	}
	
	return(prod);
}

/**********************************************************************/

void CalcM3(int n, double *beta, double **m3)
/*
  ... (comments omitted) ...
*/
{
	double *pc;
	double fr;	   /* floating point register for temp storage */
	int i;		
	double apb;	       /* sum of shape parameters of beta distribution */
	int j;	
	int nitems;
	
	apb = beta[0] + beta[1];
	nitems = n-1;
	
  /* compute first elements of all m3[j] */

	*(*(m3+n-1)) = 1.0;
	for (j=1; j<n; j++)
	{
		fr = static_cast<double>(j);
		*(*(m3+n-j-1)) = (fr / (apb+fr-1.0)) * *(*(m3+n-j));
	}

  /* compute other elements */

	for (j=1; j<n; j++)
	{
		pc = *(m3+n-j-1) + 1;
		for (i=1; i<=j; i++)
		{
			fr = apb-i+j;
			fr /= static_cast<double>(nitems - i + 1);
			fr *= *(pc-1);
			*pc++ = fr;
		}
	}
}

/*********************************************************************/

short CalcM24(int n, double para, double *m)
/*
  ... (comments omitted) ...
*/
{
	double *pm;
	double fr;	   /* floating point register for temp storage */
	int i;	
	int nitems;
	
	nitems = n-1;
	pm = m;
	*pm++ = 1.0;
	for (i=1; i<n; i++)
	{
		fr = (static_cast<double>(nitems  - (i - 1))) / (static_cast<double>(i));
		fr *= para;
		fr *= *(pm-1);
		if (fr == HUGE_VAL) return -2;	         /* check for overflow */
		*pm++ = fr;
	}

	return 0;
}

/***********************************************************************/

void CalcM15(int n, double para, double *m)
/*
  ... (comments omitted) ...
*/
{
	double *pm;
	double di;
	int i;	
	
	pm = m + (n-1);
	
	*pm-- = 1.0;
	for (i=1; i<n; i++)
	{
		di = static_cast<double>(i);
		*pm = ((para + (di - 1.0)) / di) * *(pm+1);
		--pm;
	}	
}

/*********************************************************************/

double LRChiSqr(int n, const double *rawc,
                const double *fitc, int *ncat)
/*
  ... (comments omitted) ...
*/
{
	double t;
	double chisq;
	
	chisq = 0.0;
	*ncat = 0;
	while(n--)
	{
		if (*rawc > 0.0 && *fitc > 0.0)
		{
			t = std::log(*rawc / *fitc);
			t *= *rawc;
			chisq += t;
		}
		++(*ncat);
		++fitc; ++rawc;
	}
	
	chisq *= 2.0;
	
	return(chisq);
}

/*********************************************************************/

double PChiSqr(int n, double minexpcount, const double *rawc,
	             const double *fitc, int *ncat)
/*
  ... (comments omitted) ...
*/
{
	double t;
	double chisq;
	
	double deleteraw,	             /* raw counts for categories deleted 
                                                  due fit being too small */
		     deletefit;	          /* fitted counts for deleted categories
                                               corresponding to deleteraw */
	int    deleten;		           /* number of categories deleted due to 
                                                      low expected counts */
	chisq = 0.0;
	*ncat = 0;
	deleteraw = 0.0;
	deletefit = 0.0;
	deleten = 0;
	while(n--)
	{
		if (*fitc > minexpcount)
		{
			t = *rawc - *fitc;
			t *= t;
			t /= *fitc;
			chisq += t;
			++(*ncat);
		}
		else
		{
			deleteraw += *rawc;
			deletefit += *fitc;
			++deleten;
		}
		++fitc; ++rawc;
	}
	
	/* add in chi-square value for deleted categories */

	if (deleten)
	{
		t = deleteraw - deletefit;
		t *= t;
		t /= deletefit;
		chisq += t;
		++(*ncat);
	}

	return(chisq);	
}

/*********************************************************/

double CalcLordk(double kr20, int nitems, double *rmoment)
/*
  ... (comments omitted) ...
*/
{
	double varp,dn,varr,mnm,
         k;                                 /* Lord's k */

	dn = static_cast<double>(nitems);
	
	mnm = rmoment[0] * (dn - rmoment[0]);
	varr = rmoment[1]*rmoment[1];

  /* calculate variance of item difficulties */

	varp = mnm / (dn*dn);
	varp -= (varr/dn) * (1.0 - ((dn - 1.0)/dn) * kr20);

	/* calculate k */

	k = mnm - varr - dn*varp;
	k *= 2.0;
	k = (dn*dn*(dn-1.0)*varp) / k;

  return k;
}

/*********************************************************/
