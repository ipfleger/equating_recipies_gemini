/* LogLinear.cpp  code for log-linear equating

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

  File contains functions used for:
      (a) univariate log-linear fitting under the multinomial model,
          which populates struct ULL_SMOOTH 
      (b) bivariate log-linear fitting under the multinomial model,
          which populates struct BLL_SMOOTH 

  Code follows procedures (and usually the notation) in Holland and
  Thayer (1987), which is sometimes abbreviated H&T in comments.

  A somewhat unique feature of these functions is that the user 
  can select from among various convergence criteria, including
  criteria based on moments for the actual raw scores (e.g., 
  mean, sd, skew, kurt, etc for the scores obtained using
  score() in ERutilities.c).  
*/

#include "ERutilities.h"
#include "NRutilities.h"
#include "LogLinear.h"
#include "CG_EquiEquate.h"
#include <cfloat> 
#include <cstdio>
#include <cstring>
#include <cmath>

/************************************************************************/

void Wrapper_RL(char design, char method, char smoothing,
	struct USTATS *x, struct USTATS *y,
	struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
	struct PDATA *inall, struct ERAW_RESULTS *r)
	/*
	  Wrapper for doing equipercentile equating with RG design
		and log-linear smoothing smoothing
		
	  Assumes that equating puts raw scores for x on scale of y
	  
	  NOTE: This function is used (unaltered) for both actual equating and
			equating done in Wrapper_Bootstrap().  Distinguishing between the
			two is the purpose of the variable rep

	  Input
	  
		design = 'R'(random groups)
		method = 'E'(equipercentile)
		smoothing = 'L' (log-linear smoothing)
		*x = pointer to struct USTATS (new form)
		*y = pointer to struct USTATS (old form)
		*ullx = pointer to struct BB_SMOOTH (new form)
		*ully = pointer to struct BB_SMOOTH (old form)
		rep = replication number for bootstrap; should be set to 0
			  for actual equating;
		
	  Output
		
		struct PDATA *inall:   populates selected values of inall
		
		struct ERAW_RESULTS *r: populates

		  **eraw: equated raw scores;
				  method (rows) by raw score (columns) matrix
				  of equated scores. Here there is only one method.
				  So, memory allocated for eraw[][] is:
				  eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
								 (loc(x->max,x->min,x>-inc)]
				  because we are getting equated raw scores for x
		  **mts:  moments for equated raw scores
		  
	  NOTE: If Wrapper_RL() is called in a bootstrap loop,
			then in the calling function struct ERAW_RESULTS must
			be different from struct ERAW_RESULTS for the actual
			equating.

	  Function calls other than C or NR utilities:
		EquiEquate()
		MomentsFromFD()
													
	  R. L. Brennan

	  Date of last revision: 6/30/08
	*/
{
	/* method name --- 10 characters; right justified */
	char *names[] = { "   y-equiv" };

	inall->rep = rep;               /* should be set to 0 for actual equating */
						/* counting of replications done in Wrapper_Bootstrap() */

	/* allocation and assignments for struct PDATA inall
	   Note that for every assignment of the form inall->(var) = x->(var)
	   or inall->(var) = y->(var), values vary depending on whether x or y
	   is for actual equating or a bootstrap sample; all other values are
	   the same for the actual equating and a bootstrap sample */

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

		inall->ullx = ullx;
		inall->ully = ully;
	}

	/* allocation and assignments for r */

	if (inall->rep <= 1) {         /* no storage allocation for bootstrap reps >1 */
		r->eraw = dmatrix(0, 0, 0, loc(inall->max, inall->min, inall->inc));
		r->mts = dmatrix(0, 0, 0, 3);                    /* 0,3 is for the 4 moments */
	}

	/* Compute equating results */

	EquiEquate(y->ns, y->min, y->inc, ully->crfd, x->ns, ullx->prd, r->eraw[0]);

	/* get moments */

	MomentsFromFD(inall->min, inall->max, inall->inc, r->eraw[0], inall->fdx, r->mts[0]);

}

/********************************************************************************/

void Print_RL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
	/*
	  print RL results (RG; equipercentile; log-linear smoothing)
	  
	  Input
		fp = file pointer for output
		tt[] = user supplied text identifier
		inall =  struct PDATA
		r = struct ERAW_RESULTS

	  Function calls other than C or NR utilities: None
													
	  R. L. Brennan

	  Date of last revision: 6/30/08
	*/
	int i, j;

	std::fprintf(fp, "\n\n%s\n\n", tt);

	if (inall->rep > 0)  std::fprintf(fp, "Bootstrap relication %d\n\n", inall->rep);

	std::fprintf(fp, "Equipercentile Equating with Random Groups Design\n"
		"and Polynomial Log-Linear Smoothing");

	std::fprintf(fp, "\n\nLog-Linear Smoothing for new form %c: ", inall->x->id);
	std::fprintf(fp, "\n   number of degrees of smoothing = %2d", inall->ullx->c);
	std::fprintf(fp, "\nLog-Linear Smoothing for old form %c: ", inall->y->id);
	std::fprintf(fp, "\n   number of degrees of smoothing = %2d", inall->ully->c);

	std::fprintf(fp, "\n\nInput file for new form %c: %s\n", inall->x->id, inall->xfname);
	std::fprintf(fp, "Input file for old form %c: %s\n\n", inall->y->id, inall->yfname);

	for (i = 1; i <= 30; i++) std::fprintf(fp, "-");

	/* following code set up for any number of methods */

	std::fprintf(fp, "\n\n");
	for (j = 1; j <= 15 + (inall->nm * 12 - 18) / 2; j++) std::fprintf(fp, " ");
	std::fprintf(fp, "Equated Raw Scores");
	std::fprintf(fp, "\n\n               ");
	for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "   Method %d:", j);
	std::fprintf(fp, "\nRaw Score (%c)  ", inall->x->id);
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

/*******************************************************************************/

void Wrapper_SL(char design, char method, char smoothing, struct BSTATS *xy,
	struct BLL_SMOOTH *bllxy, int rep, struct PDATA *inall,
	struct ERAW_RESULTS *r)
	/*
	  Wrapper for doing equipercentile equating with SG design
	  and log-linear smoothing.
	  
	  NOTE: This is for the SG design in which x and y do not share any items in
	  common, which means that functionally this is the external anchor case.
	  The bivariate log-linear smoothing procedure needs to know this. So, when
	  Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
	  anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
	  (with anchor set to 0) and Wrapper_SL() can still be used, but convergence
	  of the smoothing algorithm may be compromised because of dependencies between
	  x and y. (To date my experience does not suggest this is a problem.)
		
	  Assumes that equating puts raw scores for x on scale of y
	  
	  NOTE: This function is used (unaltered) for both actual equating and
			equating done in Wrapper_Bootstrap().  Distinguishing between the
			two is the purpose of the variable rep

	  Input
	  
		design = 'S' (single group)
		method = 'E' (equipercentile)
		smoothing = 'L' (log-linear smoothing)
		xy = struct BSTATS
		bllxy = struct BLL_SMOOTH
		rep = replication number for bootstrap; should be set to 0
			  for actual equating;
		
		NOTE: it is assumed that the first variable
			  in xy is indeed x and the second variable is y.
			  For bllxy, data memebers with '_x' are for x,
			  and data members with '_v' are for y.  This somewhat
			  inconsistent notation arises because BLL_SMOOTH is
			  usually used with the CG design in which the variables
			  are more naturally designated x (or u) and v.

	  Output
		
		struct PDATA inall   populates selected values of inall

		struct ERAW_RESULTS *r: populates

		  **eraw: equated raw scores;
				  method (rows) by raw score (columns) matrix
				  of equated scores. Here there is only one method.
				  So, memory allocated for eraw[][] is:
				  eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
								 (loc(x->max,x->min,x>-inc)]
				  because we are getting equated raw scores for x
		  **mts:  moments for equated raw scores
		  
	  NOTE: If Wrapper_SL() is called in a bootstrap loop,
			then in the calling function struct ERAW_RESULTS must
			be different from struct ERAW_RESULTS for the actual
			equating.
												
	  Function calls other than C or NR utilities:
		EquiEquate()
		MomentsFromFD()
													
	  R. L. Brennan

	  Date of last revision: 6/30/08
	*/
{
	/* method names --- 10 characters; right justified */
	char *names[] = { "   y-equiv" };

	inall->rep = rep;                /* should be set to 0 for actual equating */
					   /* counting of replications done in Wrapper_Bootstrap() */

	/* Allocation and assignments for struct PDATA inall>
	   Note that for every assignment of the form inall->(var) = x->(var)
	   or inall->(var) = y->(var), values vary depending on whether x or y
	   is for actual equating or a bootstrap sample; all other values are
	   the same for the actual equating and a bootstrap sample */

	if (inall->rep == 0) {     /* no assignment or stor alloc for bootstrap reps */
		std::strcpy(inall->xyfname, xy->fname);
		inall->xy = xy;
		inall->bllxy = bllxy;
		inall->design = design;
		inall->method = method;
		inall->smoothing = smoothing;
		inall->anchor = 0;  /* implicitly, anchor is external for biv log-linear
										smoothing with the SG design */

		inall->nm = 1;
		inall->names = cmatrix(0, inall->nm - 1, 0, 11);    /* only one row/method, 0 */
		std::strcpy(inall->names[0], names[0]);

		inall->min = xy->min1;
		inall->max = xy->max1;
		inall->inc = xy->inc1;
		inall->fdx = xy->fd1;
		inall->n = xy->n;
	}

	/* allocation and assignments for r */

	if (inall->rep <= 1) {         /* no storage allocation for bootstrap reps >1 */
		r->eraw = dmatrix(0, 0, 0, loc(inall->max, inall->min, inall->inc));
		r->mts = dmatrix(0, 0, 0, 3);                    /* 0,3 is for the 4 moments */
	}

	/* Compute equating results. Put x on scale of y.
	   Note that in struct xy, '1' designates x and '2' designates y;
	   in struct bllxy, '_x' designates x and '_v' designates y. So:
		 xy->ns2 = number of score categories for y
		 xy->min2 = minimum score for y
		 xy->inc2 = increment for y
		 bllxy->crfd_v = log-linear smoothed cum rel fd for y
		 xy->ns1 = number of score categories for x
		 bllxy->prd_x = log-linear smoothed PR dist for x
		 r->eraw[0] = y equivalents for x (output) */

	EquiEquate(xy->ns2, xy->min2, xy->inc2, bllxy->crfd_v, xy->ns1, bllxy->prd_x, r->eraw[0]);

	/* get moments */

	MomentsFromFD(inall->min, inall->max, inall->inc, r->eraw[0], inall->fdx, r->mts[0]);

	return;
}

/*******************************************************************************/

void Print_SL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print SL results (SG; equipercentile; log-linear smoothing)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: None
                                                
  R. L. Brennan 

  Date of last revision: 6/30/08   
*/
{
  int i,j;
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  std::fprintf(fp,"Input filename:  %s\n\n",inall->xy->fname);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
    
  std::fprintf(fp,"Equipercentile Equating with Single Group Design\n"
	  	     "and Polynomial Log-Linear Smoothing\n\n");

  std::fprintf(fp,"Bivariate Log-Linear Smoothing using the Multinomial Model\n\n");

  std::fprintf(fp,"Number of score categories for %c  = %d\n",
	         inall->xy->id1,inall->bllxy->nsu);
  std::fprintf(fp,"Number of score categories for %c  = %d\n",
	         inall->xy->id2,inall->bllxy->nsv);
  std::fprintf(fp," Total number of score categories = %d\n\n",inall->bllxy->ns);

  std::fprintf(fp,"Number of persons (i.e., total of frequencies) = %7d\n\n",
	         inall->bllxy->num_persons);

  std::fprintf(fp,"Polynomial degree for %c = %d\n",inall->xy->id1,inall->bllxy->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",inall->xy->id2,inall->bllxy->cv);
  std::fprintf(fp,"Number of cross-products = %d\n",inall->bllxy->cuv);
  if(inall->bllxy->cuv>0){
    std::fprintf(fp,"Cross-Product moments: ");
    for(i=0;i<inall->bllxy->cuv;i++) std::fprintf(fp,"  (%c%d,%c%d)%c",
                            inall->xy->id1,inall->bllxy->cpm[i][0],
							inall->xy->id2,inall->bllxy->cpm[i][1],
                            (inall->bllxy->cuv==i+1) ? ' ' : ',');
  }

  std::fprintf(fp,"\n\n");
  for(i=1;i<=30;i++) std::fprintf(fp,"-");

 /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (%c)  ",inall->xy->id1);
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"  %s",inall->names[j]);
  std::fprintf(fp,"\n\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    std::fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  std::fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][0]);
  std::fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][1]);
  std::fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][2]);
  std::fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][3]);
  
  std::fprintf(fp,"\n\n");
  for(j=1;j<=63;j++) std::fprintf(fp,"*");

  return;

}

/*******************************************************************************/

void Wrapper_CL(char design, char method, char smoothing, 
                double w1, int anchor, double rv1, double rv2,
                struct BSTATS *xv, struct BSTATS *yv, 
				struct BLL_SMOOTH *bllxv, struct BLL_SMOOTH *bllyv,
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for equipercentile equating for CG design with log-linear smoothing. 
  Equipercentile equating includes frequency estimation with 
  Braun-Holland (linear) results, modified frequency estimation with 
  Braun-Holland (linear) results, and chained equipercentile equating
    
  Assumes that in xv, score 1 is for x and score 2 is for v
  Assumes that in yv, score 1 is for y and score 2 is for v
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep
  
  Input:
  
    design = 'C' (CINEG)

    method:  'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
             'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
             'G' = FE + BH-FE + MFE + BH-MFE
             'C' = Chained
			 'H' = FE + BH-FE + Chained
             'A' = FE + BH-FE + MFE + BH-MFE + Chained
              
    smoothing = 'L' (log-linear) 

    w1 = weight for pop. 1 (associated with xv)
         [0,1] except that for any number outside this 
         range, proportional weights are used -- i.e.,
         w1 = xv->n/(xv->n + yv->n)
    anchor = 0 --> external; otherwise internal
    rv1 = reliability of common items for population 1 
          (set to 0 for all methods except 'F', 'G, and 'A')
    rv2 = reliability of common items for population 2
          (set to 0 for all methods except 'F', 'G, and 'A')
    xv = struct BSTATS
    yv = struct BSTATS 
    bllxv = struct BLL_SMOOTH; uses brfd, prd_x,  and crfd_v 
	bllyv = struct BLL_SMOOTH; uses brfd, crfd_x, and crfd_v
	        (Note that bllyv->crfd_x is really crfd for y in pop 2)
    rep = replication number for bootstrap; should be set to 0
          for actual equating; 
   
    NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
    
  Output:
    
    struct PDATA inall:   populates selected values of inall 
    
    struct ERAW_RESULTS r          
 
      a[] = slopes for Braun-Holland
      b[] = intercepts for Braun-Holland
      eraw[][]:  equated raw scores
      mts[][]:  moments for equated raw scores   
          
      NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
            of equated scores; memory allocated here;
            eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
            because we are getting equated raw scores for x.
            eraw[][] stored in this "row"  manner so that 
            Equated_ss() can be used directly on 
            eraw[k] where k is the method number  
          
  NOTE: Whenever method differs, there must be different structures
        passed as struct PDATA and struct ERAW_RESULTS 
    
  NOTE: If Wrapper_CL() is called in a bootstrap loop, then in
        the calling function struct ERAW_RESULTS must be different
        from struct ERAW_RESULTS for the actual equating. 
                                            
  Function calls other than C or NR utilities:
    FEorMFE_EE()
    Chained_EE()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08   
*/
{ 
  int i;
  double *ptr;                                      /* pointer for eraw[] */
                       /* method names --- 10 characters; right justified */
  char *names[] ={"        FE", "     BH-FE", "       MFE", "    BH-MFE",
                  "  ChainedE"};                   
  
  inall->rep = rep;               /* should be set to 0 for actual equating. */
                    /* Counting of replications done in Wrapper_Bootstrap(), 
             which is why this statement cannot be in the if statement below */ 
                    
  /* allocation and assignments for inall
     Note that for every assignment of the form inall->(var) = xv->(var)
	 or inall->(var) = yv->(var) values vary depending on whether xv or yv
	 is for actual equating or a bootstrap sample; all other values are 
	 the same for the actual equating and a bootstrap sample */
  
  if(inall->rep == 0){   /* no assignment or stor alloc for bootstrap reps */
    std::strcpy(inall->xfname,xv->fname);
    std::strcpy(inall->yfname,yv->fname);
    inall->xv = xv;
    inall->yv = yv;
    inall->bllxv = bllxv;
	inall->bllyv = bllyv;
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
    inall->w1 = (w1<0 || w1>1) ? static_cast<double>(xv->n)/(xv->n + yv->n) : w1; 
                                   /* proportional wts if w1 outside [0,1] */          
    inall->anchor = anchor;
    inall->rv1 = rv1;
    inall->rv2 = rv2;

    if((method == 'F' || method == 'G' || method == 'A') && 
       (rv1 == 0 || rv2 == 0))
       runerror("\nMFE cannot be conducted since rv1 == 0 or rv2 == 0");
    
    inall->names = cmatrix(0,4,0,11);             /* maximum of five names */
 
    if(method == 'E'){                                    /* method == 'E' */
      inall->nm = 2;
      for(i=0;i<=1;i++) std::strcpy(inall->names[i],names[i]);
    }
    else if(method == 'F'){                               /* method == 'F' */
      inall->nm = 2;
      for(i=2;i<=3;i++) std::strcpy(inall->names[i-2],names[i]);
    }
    else if(method == 'G'){                               /* method == 'G' */
      inall->nm = 4;
      for(i=0;i<=3;i++) std::strcpy(inall->names[i],names[i]);
    }
    else if(method == 'C'){                               /* method == 'C' */
      inall->nm = 1;
      std::strcpy(inall->names[0],names[4]);
    }
	else if(method == 'H'){                               /* method == 'H' */
      inall->nm = 3;
	  for(i=0;i<=1;i++) std::strcpy(inall->names[i],names[i]);
      std::strcpy(inall->names[2],names[4]);
    }
    else{                                                 /* method == 'A' */
      inall->nm = 5;
      for(i=0;i<=4;i++) std::strcpy(inall->names[i],names[i]);
    }  

    inall->min = xv->min1;  
    inall->max = xv->max1;
    inall->inc = xv->inc1;
    inall->fdx = xv->fd1;
    inall->n = xv->n;
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){         /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,inall->nm-1,0,loc(xv->max1,xv->min1,xv->inc1)); 
    r->mts = dmatrix(0,inall->nm-1,0,3);          /* 0,3 is for the 4 moments */
    r->fxs = dmatrix(0,1,0,loc(xv->max1,xv->min1,xv->inc1));
    r->gys = dmatrix(0,1,0,loc(yv->max1,yv->min1,yv->inc1));
  }
   
/* Equipercentile results, including Braun-Holland (BH) linear. 
   Note: For FE syn densities are in fxs[0] and gys[0]
         For MFE syn densities are in fxs[1] and gys[1] 
         For BH under FE, slope in a[0] and intercept in b[0]
         For BH under MFE, slope in a[1] and intercept in b[1] */

                                           /* FE + BH-FE in positions 0 and 1*/
  if(method == 'E' || method == 'G' || method == 'A' || method == 'H')     
    FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1, 
               yv->ns1, yv->min1, yv->max1, yv->inc1,
               bllxv->brfd, bllyv->brfd, 0, 0, 
               r->fxs[0], r->gys[0], r->eraw[0],
               &r->a[0], &r->b[0], r->eraw[1]); 

  if(method == 'F')                     /* MFE + BH-MFE in positions 0 and 1 */
    FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1, 
               yv->ns1, yv->min1, yv->max1, yv->inc1,
               bllxv->brfd, bllyv->brfd, inall->rv1, inall->rv2, 
               r->fxs[1], r->gys[1], r->eraw[0],
               &r->a[1], &r->b[1], r->eraw[1]);    

  if(method == 'G' || method == 'A')    /* MFE + BH-MFE in positions 2 and 3 */
    FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1, 
               yv->ns1, yv->min1, yv->max1, yv->inc1,
               bllxv->brfd, bllyv->brfd, inall->rv1, inall->rv2, 
               r->fxs[1], r->gys[1], r->eraw[2],
               &r->a[1], &r->b[1], r->eraw[3]); 

  /* Chained equipercentile method. Note that smoothing is bivariate
     log-linear smoothing, not univariate log-linear smoothing.
	 if(method == 'C')  Chained in position 0 
	 if(method == 'A')  Chained in position 4 
	 if(method == 'H')  Chained in position 2 */

  if(method == 'C') ptr = r->eraw[0];
  else if(method == 'A') ptr = r->eraw[4];
  else if(method == 'H') ptr = r->eraw[2];
  else ptr = nullptr;

  if(ptr)
    Chained_EE(xv->ns1, 
	             bllxv->prd_x,    /* this is smoothed pr dist for x in pop 1 */ 
			   xv->min2, 
			   xv->max2, 
			   xv->inc2, 
               xv->ns2, 
			     bllxv->crfd_v,      /* this is smoothed crfd for v in pop 1 */
			   yv->min1, 
			   yv->inc1, 
			   yv->ns1,   
                 bllyv->crfd_x,      /* this is smoothed crfd for y in pop 2 */ 
			     bllyv->crfd_v,      /* this is smoothed crfd for v in pop 2 */
				 ptr);  
                        
/* get moments */

  for(i=0;i<=inall->nm-1;i++) 
    MomentsFromFD(xv->min1,xv->max1,xv->inc1,r->eraw[i],xv->fd1,r->mts[i]);
  
  return;
}

/****************************************************************************/

void Print_CL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print CL results (CG; log-linear smoothed equipercentile equating)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: 
    score() 
                                            
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
{
  int i,j,
      jlast = (inall->method=='A') ? 75 : 63;

  char x = inall->xv->id1,                            /* x = u[] + v for new form */
       v = inall->xv->id2,                           /* common items for new form */
       y = inall->yv->id1,                            /* y = w[] + z for old form */
       z = inall->yv->id2,                           /* common items for old form */
       u[3],                  /* string designating non-common items for new form */
       w[3];                  /* string designating non-common items for old form */

  if(inall->anchor){u[0] = x; u[1] = '\''; u[2] = '\0';}  /* internal anchor "x'" */
  else {u[0] = ' '; u[1] = x; u[2] = '\0';}               /* external anchor " x" */

  if(inall->anchor){w[0] = y; w[1] = '\''; w[2] = '\0';}  /* internal anchor "y'" */
  else {w[0] = ' '; w[1] = y; w[2] = '\0';}               /* external anchor " y" */
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  std::fprintf(fp,"Equipercentile Equating with CINEG Design"
             " (%s Anchor; w1 = %7.5f)\n"
			 "and Polynomial Log-Linear Smoothing\n\n",
             (inall->anchor != 0) ? "Internal" : "External", inall->w1); 
  std::fprintf(fp,"Input file for %c%c and pop 1: %s\n",x,v,inall->xfname);
  std::fprintf(fp,"Input file for %c%c and pop 2: %s\n\n",y,z,inall->yfname);

  if(inall->anchor)
	  std::fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score; and\n"
				 "(2) variables for old form and pop 2 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score.\n\n",
			     x,u,v,x,u,v,x,u,v,  y,w,z,y,w,z,y,w,z);
  else
	  std::fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
				 "    %c and %c, where %c does not include %c;\n"
				 "(2) variables for old form and pop 2 are identified as\n"
				 "    %c and %c, where %c does not include %c.\n\n",
			     x,v,x,v, y,z,y,z);

  /* information about x and v in pop 1 */

  std::fprintf(fp,"Number of persons for %c = %7d\n\n",x,inall->bllxv->num_persons);
  std::fprintf(fp,"Number of score categories for %c = %d\n",x,inall->bllxv->nsx);
  std::fprintf(fp,"Number of score categories for %c = %d\n",v,inall->bllxv->nsv);

  std::fprintf(fp,"Polynomial degree for %s = %d\n",u,inall->bllxv->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",v,inall->bllxv->cv);
  std::fprintf(fp,"Number of cross-products %s%c = %d\n",u,v,inall->bllxv->cuv);
  if(inall->bllxv->cuv>0) std::fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllxv->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               u,inall->bllxv->cpm[i][0],v,inall->bllxv->cpm[i][1],
                               (inall->bllxv->cuv==i+1) ? ' ' : ',');

 /* information about y and v in pop 2 */

  std::fprintf(fp,"\n\nNumber of persons for %c = %7d\n\n",y,inall->bllyv->num_persons);
  std::fprintf(fp,"Number of score categories for %c = %d\n",y,inall->bllyv->nsx);
  std::fprintf(fp,"Number of score categories for %c = %d\n",z,inall->bllyv->nsv);

  std::fprintf(fp,"Polynomial degree for %s = %d\n",w,inall->bllyv->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",z,inall->bllyv->cv);
  std::fprintf(fp,"Number of cross-products %s%c = %d\n",w,z,inall->bllyv->cuv);
  if(inall->bllyv->cuv>0) std::fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllyv->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               w,inall->bllyv->cpm[i][0],z,inall->bllyv->cpm[i][1],
                               (inall->bllyv->cuv==i+1) ? ' ' : ',');
  
  /* print equated raw scores for methods */

  if(inall->method == 'E' || inall->method == 'G' || inall->method == 'A')
    std::fprintf(fp,"\n\nBraun-Holland (BH) Under Frequency Estimation (FE):\n"
               "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
               r->b[0], r->a[0]);
  if(inall->method == 'F' || inall->method == 'G' || inall->method == 'A')
    std::fprintf(fp,"Braun-Holland (BH) Under Modified Frequency Estimation (MFE):\n"
               "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
               r->b[1], r->a[1]);

  if(inall->method == 'E' || inall->method == 'F' || 
	 inall->method == 'G' || inall->method == 'A')
     for(i=1;i<=jlast;i++) std::fprintf(fp,"-");

  /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (%c)  ",x);
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"  %s",inall->names[j]);
  std::fprintf(fp,"\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    std::fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  std::fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][0]);
  std::fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][1]);
  std::fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][2]);
  std::fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][3]);
  
  std::fprintf(fp,"\n\n");
  for(j=1;j<=jlast;j++) std::fprintf(fp,"*");

}

/****************************************************************************/

void Wrapper_Smooth_ULL(struct USTATS *x, int c,
						int scale, int Btype, int ctype, double crit,
						FILE *fp, struct ULL_SMOOTH *s)
/*
  Wrapper to do univariate log-linear smoothing.

  Input
    x   =  UTATS structure
    c   = number of degrees for polynomial smoothing
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    fp  = pointer to output file (if NULL then no output printed;
	      in particular, results for each iteration step are not printed)

  Output
    populates s, which is a ULL_SMOOTH structure

  Function calls other than C or NR utilities:  
    Smooth_ULL()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  Smooth_ULL(x->n, x->ns, x->min, x->inc, x->dbl_fd, c, 
	         scale, Btype, ctype, crit, fp, s); 

  return; 
} 

/**************************************************************/

void Smooth_ULL(int n, int ns, double min, double inc,
                double *fd, int c, 
				int scale, int Btype, int ctype, double crit,
				FILE *fp, struct ULL_SMOOTH *s)
/*
  Performs univariate log-linear smoothing in terms of the 
  multinomial model as described by Holland and Thayer (1987),
  abbreviated here as H&T.  

  Input
  
    n = number of persons
    ns = number of score categories
    min = minimum raw score
    inc = increment in raw scores
    fd[] = frequency distribution
    c = number of degrees for polynomial smoothing
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    *fp = pointer to output file
    
  Output: populates struct ULL_SMOOTH s 

  Function calls other than C or NR utilities:
    design_matrix()
    iteration()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,
      max_nit = 40;           /* maximum number of iterations */ 

  s->num_persons = n;
  s->ns = ns;
  s->min = min;
  s->inc = inc;
  s->c = c;
  s->nct = fd;
  s->scale = scale;
  s->max_nit = max_nit;
  s->ctype = ctype;
  s->Btype = Btype;
  s->crit = crit;

  s->B_raw = dmatrix(0,ns-1,0,c-1);
  s->B = dmatrix(0,ns-1,0,c-1);
  s->mct = dvector(0,ns-1);
  s->Beta = dvector(0,c-1);
  s->n_mts = dvector(0,c-1);
  s->m_mts = dvector(0,c-1);
  s->n_mts_raw = dvector(0,c-1);
  s->m_mts_raw = dvector(0,c-1);

  s->density = dvector(0,ns-1);
  s->crfd = dvector(0,ns-1);
  s->prd = dvector(0,ns-1);
                      
  design_matrix(ns,min,inc, 0,0,0, c,0,0,nullptr, scale,s->B_raw,s->B); 

  /* iteration-step results not printed if first parameter is NULL;
     results are printed if first parameter is fp!=NULL. Note that fp
	 can be set to NULL in Wrapper_Smooth_ULL() */ 

  s->nit = iteration(fp, s->B, s->B_raw, fd, static_cast<double>(n), 
                     nullptr, ns, c, 0, 0, nullptr,
                     max_nit, ctype, Btype, crit,
                     s->Beta, s->mct, s->n_mts, s->m_mts, 
                     s->n_mts_raw, s->m_mts_raw, 
                     &(s->lrchisq),&(s->nzero),&(s->ap));

  for (i=0;i<ns;i++) s->density[i] = s->mct[i]/s->num_persons;
  cum_rel_freqs(0, ns-1, 1, s->density, s->crfd); 
  for (i=0;i<ns;i++)
    s->prd[i] = perc_rank(0, ns-1, 1, s->crfd, static_cast<double>(i));

  return;
}
/****************************************************************************/

void Print_CL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print CL results (CG; log-linear smoothed equipercentile equating)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: 
    score() 
                                            
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
{
  int i,j,
      jlast = (inall->method=='A') ? 75 : 63;

  char x = inall->xv->id1,                            /* x = u[] + v for new form */
       v = inall->xv->id2,                           /* common items for new form */
       y = inall->yv->id1,                            /* y = w[] + z for old form */
       z = inall->yv->id2,                           /* common items for old form */
       u[3],                  /* string designating non-common items for new form */
       w[3];                  /* string designating non-common items for old form */

  if(inall->anchor){u[0] = x; u[1] = '\''; u[2] = '\0';}  /* internal anchor "x'" */
  else {u[0] = ' '; u[1] = x; u[2] = '\0';}               /* external anchor " x" */

  if(inall->anchor){w[0] = y; w[1] = '\''; w[2] = '\0';}  /* internal anchor "y'" */
  else {w[0] = ' '; w[1] = y; w[2] = '\0';}               /* external anchor " y" */
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  std::fprintf(fp,"Equipercentile Equating with CINEG Design"
             " (%s Anchor; w1 = %7.5f)\n"
			 "and Polynomial Log-Linear Smoothing\n\n",
             (inall->anchor != 0) ? "Internal" : "External", inall->w1); 
  std::fprintf(fp,"Input file for %c%c and pop 1: %s\n",x,v,inall->xfname);
  std::fprintf(fp,"Input file for %c%c and pop 2: %s\n\n",y,z,inall->yfname);

  if(inall->anchor)
	  std::fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score; and\n"
				 "(2) variables for old form and pop 2 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score.\n\n",
			     x,u,v,x,u,v,x,u,v,  y,w,z,y,w,z,y,w,z);
  else
	  std::fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
				 "    %c and %c, where %c does not include %c;\n"
				 "(2) variables for old form and pop 2 are identified as\n"
				 "    %c and %c, where %c does not include %c.\n\n",
			     x,v,x,v, y,z,y,z);

  /* information about x and v in pop 1 */

  std::fprintf(fp,"Number of persons for %c = %7d\n\n",x,inall->bllxv->num_persons);
  std::fprintf(fp,"Number of score categories for %c = %d\n",x,inall->bllxv->nsx);
  std::fprintf(fp,"Number of score categories for %c = %d\n",v,inall->bllxv->nsv);

  std::fprintf(fp,"Polynomial degree for %s = %d\n",u,inall->bllxv->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",v,inall->bllxv->cv);
  std::fprintf(fp,"Number of cross-products %s%c = %d\n",u,v,inall->bllxv->cuv);
  if(inall->bllxv->cuv>0) std::fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllxv->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               u,inall->bllxv->cpm[i][0],v,inall->bllxv->cpm[i][1],
                               (inall->bllxv->cuv==i+1) ? ' ' : ',');

 /* information about y and v in pop 2 */

  std::fprintf(fp,"\n\nNumber of persons for %c = %7d\n\n",y,inall->bllyv->num_persons);
  std::fprintf(fp,"Number of score categories for %c = %d\n",y,inall->bllyv->nsx);
  std::fprintf(fp,"Number of score categories for %c = %d\n",z,inall->bllyv->nsv);

  std::fprintf(fp,"Polynomial degree for %s = %d\n",w,inall->bllyv->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",z,inall->bllyv->cv);
  std::fprintf(fp,"Number of cross-products %s%c = %d\n",w,z,inall->bllyv->cuv);
  if(inall->bllyv->cuv>0) std::fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllyv->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               w,inall->bllyv->cpm[i][0],z,inall->bllyv->cpm[i][1],
                               (inall->bllyv->cuv==i+1) ? ' ' : ',');
  
  /* print equated raw scores for methods */

  if(inall->method == 'E' || inall->method == 'G' || inall->method == 'A')
    std::fprintf(fp,"\n\nBraun-Holland (BH) Under Frequency Estimation (FE):\n"
               "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
               r->b[0], r->a[0]);
  if(inall->method == 'F' || inall->method == 'G' || inall->method == 'A')
    std::fprintf(fp,"Braun-Holland (BH) Under Modified Frequency Estimation (MFE):\n"
               "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
               r->b[1], r->a[1]);

  if(inall->method == 'E' || inall->method == 'F' || 
	 inall->method == 'G' || inall->method == 'A')
     for(i=1;i<=jlast;i++) std::fprintf(fp,"-");

  /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (%c)  ",x);
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"  %s",inall->names[j]);
  std::fprintf(fp,"\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    std::fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  std::fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][0]);
  std::fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][1]);
  std::fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][2]);
  std::fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][3]);
  
  std::fprintf(fp,"\n\n");
  for(j=1;j<=jlast;j++) std::fprintf(fp,"*");

}

/****************************************************************************/

void Wrapper_Smooth_ULL(struct USTATS *x, int c,
						int scale, int Btype, int ctype, double crit,
						FILE *fp, struct ULL_SMOOTH *s)
/*
  Wrapper to do univariate log-linear smoothing.

  Input
    x   =  UTATS structure
    c   = number of degrees for polynomial smoothing
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    fp  = pointer to output file (if NULL then no output printed;
	      in particular, results for each iteration step are not printed)

  Output
    populates s, which is a ULL_SMOOTH structure

  Function calls other than C or NR utilities:  
    Smooth_ULL()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  Smooth_ULL(x->n, x->ns, x->min, x->inc, x->dbl_fd, c, 
	         scale, Btype, ctype, crit, fp, s); 

  return; 
} 

/**************************************************************/

void Smooth_ULL(int n, int ns, double min, double inc,
                double *fd, int c, 
				int scale, int Btype, int ctype, double crit,
				FILE *fp, struct ULL_SMOOTH *s)
/*
  Performs univariate log-linear smoothing in terms of the 
  multinomial model as described by Holland and Thayer (1987),
  abbreviated here as H&T.  

  Input
  
    n = number of persons
    ns = number of score categories
    min = minimum raw score
    inc = increment in raw scores
    fd[] = frequency distribution
    c = number of degrees for polynomial smoothing
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    *fp = pointer to output file
    
  Output: populates struct ULL_SMOOTH s 

  Function calls other than C or NR utilities:
    design_matrix()
    iteration()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,
      max_nit = 40;           /* maximum number of iterations */ 

  s->num_persons = n;
  s->ns = ns;
  s->min = min;
  s->inc = inc;
  s->c = c;
  s->nct = fd;
  s->scale = scale;
  s->max_nit = max_nit;
  s->ctype = ctype;
  s->Btype = Btype;
  s->crit = crit;

  s->B_raw = dmatrix(0,ns-1,0,c-1);
  s->B = dmatrix(0,ns-1,0,c-1);
  s->mct = dvector(0,ns-1);
  s->Beta = dvector(0,c-1);
  s->n_mts = dvector(0,c-1);
  s->m_mts = dvector(0,c-1);
  s->n_mts_raw = dvector(0,c-1);
  s->m_mts_raw = dvector(0,c-1);

  s->density = dvector(0,ns-1);
  s->crfd = dvector(0,ns-1);
  s->prd = dvector(0,ns-1);
                      
  design_matrix(ns,min,inc, 0,0,0, c,0,0,nullptr, scale,s->B_raw,s->B); 

  /* iteration-step results not printed if first parameter is NULL;
     results are printed if first parameter is fp!=NULL. Note that fp
	 can be set to NULL in Wrapper_Smooth_ULL() */ 

  s->nit = iteration(fp, s->B, s->B_raw, fd, static_cast<double>(n), 
                     nullptr, ns, c, 0, 0, nullptr,
                     max_nit, ctype, Btype, crit,
                     s->Beta, s->mct, s->n_mts, s->m_mts, 
                     s->n_mts_raw, s->m_mts_raw, 
                     &(s->lrchisq),&(s->nzero),&(s->ap));

  for (i=0;i<ns;i++) s->density[i] = s->mct[i]/s->num_persons;
  cum_rel_freqs(0, ns-1, 1, s->density, s->crfd); 
  for (i=0;i<ns;i++)
    s->prd[i] = perc_rank(0, ns-1, 1, s->crfd, static_cast<double>(i));

  return;
}

/*******************************************************************/

void Print_ULL(FILE *fp, char tt[], struct USTATS *x,
               struct ULL_SMOOTH *s, int print_dm, int print_mts)
/*
  print results in ULL_SMOOTH s
  
  Input
    fp   = file pointer for output
    tt[] = user supplied text identifier
    x    = struct USTATS 
    s    = struct ULL_SMOOTH
    print_dm =  print design matrices
	            0 --> no printing of matrices 
	            1 --> print B (and coefficients of B)
                2 --> print B and B_raw (and coefficients of B)
				NOTE: if scale = 0, then only B_raw is printed
   print_mts = print moments
	            0 --> no print
				1 --> print moments for B
				2 --> print moments for B and B-raw
				NOTE: if scale = 0, then only moments of B_raw are printed
                   					
  Function calls other than C or NR utilities: 
    score()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j;

  std::fprintf(fp,"");

  std::fprintf(fp,"\n\n%s\n\n",tt);
  std::fprintf(fp,"Input filename:  %s\n\n",x->fname);

  std::fprintf(fp,"Univariate Log-Linear Smoothing using the Multinomial Model\n\n");

  std::fprintf(fp,"Number of score categories = %d\n",s->ns);
  std::fprintf(fp,"Number of persons (i.e., total of frequencies) = %7d\n\n",
	         s->num_persons);

  std::fprintf(fp,"Polynomial degree = %d\n\n",s->c);

  if(s->scale==0) std::fprintf(fp,"Design matrix B is powers of raw scores\n\n");
  else std::fprintf(fp,"Design matrix B is scaled powers of raw scores\n"
                  "such that for any column of B, sum(elements) = 0 and\n"
                  "sum(elements^2) = 1\n\n");
  
  std::fprintf(fp,"Convergence criterion:  ");
  if(s->Btype==0 && s->ctype==0)
    std::fprintf(fp,"|(n_mts[i] - m_mts[i]| <= %15.10f",s->crit); 
  else if(s->Btype==0 && s->ctype==1)
    std::fprintf(fp,"|n_mts[i] - m_mts[i])/n_mts[i]| <= %15.10f",s->crit);
  else if(s->Btype==1 && s->ctype==0)
    std::fprintf(fp,"|(n_mts_raw[i] - m_mts_raw[i]| <= %15.10f",s->crit); 
  else 
    std::fprintf(fp,"|n_mts_raw[i] - m_mts_raw[i])/n_mts_raw[i]|"
                 " <= %15.10f",s->crit);
  std::fprintf(fp,"  for i = 1 to %d\n\n",s->c);

  std::fprintf(fp,"Number of iterations to convergence = %d\n\n",s->nit);

  std::fprintf(fp,"Likelihood-ratio chi-square = %12.5f "
             "with %d degrees of freedom\n\n",s->lrchisq,s->ns-s->c-1);

  /* NOTE: if scale=0, the B_raw matrix is the design matric B.  In 
     this case only B_raw is printed along with its coefficients
	 and mts_raw[] */

  if(s->scale==0) std::fprintf(fp,"Since scale=0, B_raw=B is the design matrix\n\n");

  /* print coefficients for design matrix */

  std::fprintf(fp,"COEFFICIENTS:\n\n");
    for(i=1;i<=s->c;i++)
      std::fprintf(fp,"   Beta[%2d] = %12.5f\n",i,s->Beta[i-1]); 

  /* print B matrix (scaled) if print_dm is positive and scale==1 */

  if(print_dm && s->scale==1){ 

    std::fprintf(fp,"\n\nB DESIGN MATRIX");
    std::fprintf(fp,"\n\nCategory   Score");
    for(j=0;j<s->c;j++) std::fprintf(fp,"    Deg %2d",j+1); 
    std::fprintf(fp,"\n"); 
    for(i=0;i<s->ns;i++){
      std::fprintf(fp,"\n%8d%8.3f",i,score(i,s->min,s->inc));
      for(j=0;j<s->c;j++)
        std::fprintf(fp,"%10.5f",s->B[i][j]);
    }
  }

  /* print B_raw matrix if print_dm==2 and scale==1,
     or print_dm==1 and scale==0 */

  if( (print_dm==2 && s->scale==1)  || (print_dm && s->scale==0) ){

    std::fprintf(fp,"\n\nB_raw DESIGN MATRIX");
    std::fprintf(fp,"\n\nCategory   Score");
    for(j=0;j<s->c;j++) std::fprintf(fp,"           Degree %2d",j+1);
    std::fprintf(fp,"\n"); 
    for(i=0;i<s->ns;i++){
      std::fprintf(fp,"\n%8d%8.3f",i,score(i,s->min,s->inc));
      for(j=0;j<s->c;j++)
        std::fprintf(fp,"%20.5f",s->B_raw[i][j]);
    }
  }

  /* print actual and fitted frequencies */

  std::fprintf(fp,"\n\n");
  for(i=1;i<=93;i++) std::fprintf(fp,"-");

  std::fprintf(fp,"\n\n                      Actual      Fitted\n");
  std::fprintf(fp,"                      freqs:      freqs:\n");
  std::fprintf(fp,"Category   Score       nct[]       mct[]          prop  fitted-prop"
	         "  fitted-crfd   fitted-prd\n\n");
  for(i=0;i<s->ns;i++)
    std::fprintf(fp,"%8d%8.3f%12.0f%12.5f       %7.5f %12.5f %12.5f %12.5f\n",
	        i,score(i,s->min,s->inc),s->nct[i],s->mct[i],
			s->nct[i]/s->num_persons,s->density[i],
			s->crfd[i],s->prd[i]);

  /* print moments of B matrix (scaled) if print_mts is positive and scale==1 */

  if(print_mts && s->scale==1){
    std::fprintf(fp,"\n");
    for(i=1;i<=s->c;i++)
      std::fprintf(fp,"         mts[%2d]%12.5f%12.5f\n",
                i,s->n_mts[i-1], s->m_mts[i-1]);  
  }

  /* print B_raw moments if print_mts==2 and scale==1,
     or print_mts==1 and scale==0 */

  if( (print_mts==2 && s->scale==1)  || (print_mts && s->scale==0) ){
    std::fprintf(fp,"\n");
    for(i=1;i<=s->c;i++)
      std::fprintf(fp,"     mts_raw[%2d]%12.5f%12.5f\n",
                 i,s->n_mts_raw[i-1], s->m_mts_raw[i-1]);  
  }

  if(print_mts) std::fprintf(fp,"\nNOTES:\n\n");

  if(print_mts && s->scale==1){	
    std::fprintf(fp,  
      " mts[] are moments associated with the design matrix B that is\n"
      " used to obtain the log-linear solution. Columns of B are either:\n" 
      " (a) powers of raw scores (using min and inc); this is\n"  
      "     accomplished by setting scale=0 in Smooth_BLL(); or\n" 
      " (b) scaled powers of raw scores such that for each column of\n"  
      "     B, sum(elements) = 0 and sum(elements^2) = 1; this is\n"
      "     accomplished by setting scale=1 in Smooth_BLL().\n" 
      " (Experience suggests that scale=1 leads to more numerical\n"    /* edit change 3-19-09 */
      " stability.) Letting i be a score category, the jth moment is\n" /* edit change 3-19-09 */
      " mts[j] = sum_i(B[i][j-1]*rf[i]) where rf[] is either actual\n" 
      " relative frequencies or fitted relative frequencies.\n\n");
  }

   if( (print_mts==2 && s->scale==1)  || (print_mts && s->scale==0) ){
    std::fprintf(fp,
      " mts_raw[] are based on (a) above (i.e, no scaling). Importantly,\n" 
      " however, the moments greater than or equal to 3 are central moments\n" 
      " associated with matrix B_raw in the following sense. Letting i be a\n"
      " score category, rf[i] be a relative frequency, and mean = sum_i(B_raw[i][0]*rf[i]),\n"
	  " then:\n"
      " mts_raw[1] = mean,  mts_raw[2] = sd, and for j = 3,...,%d,\n" 
      " mts_raw[j] =  sum_i(B_raw[i][0] - mean)^(j)*rf[i]/sd^(j).\n\n",s->c);
  }
  
  for(i=1;i<=70;i++)  std::fprintf(fp,"*");

  return;
} 

/*******************************************************************/ 

/*******************************************************************/ 

void Wrapper_Smooth_BLL(struct BSTATS *xv, int anchor,
                        int cu, int cv, int cuv, int cpm[][2],
						int scale, int Btype, int ctype, double crit,
                        FILE *fp, struct BLL_SMOOTH *s)
/*
  Wrapper to do univariate log-linear smoothing.

  NOTE: if this function is used with the SG design,
        anchor should be set to 0 (external)

  Input
    xv      =  BTATS structure
    anchor  = 0 (external); = 1 (internal) 
    cu      = number of degrees of smoothing for u
    cv      = number of degrees of smoothing for v
    cuv     = number of cross-product moments
    cpm[cuv-1][2] = zero-offset matrix designating cross-
                    product moments.  Example: let cuv = 3,
                    and the desired cross-product moments be 
                    (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                    Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                    cpm[2] = {2,1}. 
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    fp  = pointer to output file (if NULL then no output printed;
	      in particular, results for each iteration step are not printed) 

  Output
    gets s->nct and 
    calls function that populates other elements of s

  NOTE: When anchor is internal, smoothing is done on matrix of
       (non-common-item scores) by (common-item scores).  For an
	   internal anchor, smoothing could be done on matrix of
	   (total scores) by (common-item scores), and convergence
	   should not be a problem, but structural zeros likely would be
	   given fitted positive values.  (For an internal anchor,
	   structural zeros occur when a total score is impossible
	   given a particular common item score.)

  Function calls other than C or NR utilities:  
    get_nct_bfd()
    Smooth_BLL()
    get_bfd_mct()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j,
      nsx = xv->ns1, 
      nsv = xv->ns2,
	  ns = nsx*nsv,
      nsu = (anchor) ? xv->ns1 - xv->ns2 + 1: xv->ns1;
  double *ptr,                                /* pointer variable */
	     minu = (anchor) ?  xv->min1 - xv->min2: xv->min1,
         incu = xv->inc1;    /* to be sensible, it must be true that
                                       xv->inc1 = xv->inc2 = incu */

  /* nsu is number of categories for non-common items */
 
  s->nct = dvector(0,nsu*nsv-1);
  s->bfd = dmatrix(0,nsx-1,0,nsv-1);  

  /* For an external anchor, directly convert xv->bfd[][] 
     to a row-major vector nct[]. For an internal anchor, 
     conceptually we first collapse xv->bfd[][] to 
     uv->bfd[][], and then convert to nct[].  
     In the end this means that nct[] is a row-major vector 
     with nsu*nsv elements */

  get_nct_bfd(anchor, xv->ns1, xv->ns2, xv->dbl_bfd, s->nct); 

  /* call bivariate log-linear smoothing function, which
     assumes that row categores are NOT included in
     column categories */ 

  Smooth_BLL(xv->n, nsu, minu, incu, 
             xv->ns2, xv->min2, xv->inc2, s->nct, 
             anchor, cu, cv, cuv, cpm, 
			 scale, Btype, ctype, crit, fp, s); 

  /* For an external anchor, directly convert the row major 
     vector s->mct[nsu*nsv] to s->bfd[nsu][nsv]. For an 
     internal anchor, convert s->mct[nsu*nsv] to
     s->bfd[nsu+nsv-1][nsv] repositioning elements and
	 adding structural zeros. */

  get_bfd_mct(anchor, xv->ns1, xv->ns2, s->mct, s->bfd);

  /* get row and col marginal densities, crfd's and prd's for bfd[][] */

  s->nsx = nsx;
  s->minx = xv->min1;
  s->incx = xv->inc1;

  s->fd_x = dvector(0,nsx-1);    /* row marginal frequencies for bfd[][] */
  s->density_x = dvector(0,nsx-1);   /* row marginal density for bfd[][] */
  s->crfd_x = dvector(0,nsx-1);         /* row marginal crfd for bfd[][] */
  s->prd_x = dvector(0,nsx-1);            /* row marginal PR for bfd[][] */

  s->fd_v = dvector(0,nsv-1);    /* col marginal frequencies for bfd[][] */
  s->density_v = dvector(0,nsv-1);   /* col marginal density for bfd[][] */
  s->crfd_v = dvector(0,nsv-1);         /* col marginal crfd for bfd[][] */
  s->prd_v = dvector(0,nsv-1);            /* col marginal PR for bfd[][] */
                      
  for(j=0;j<nsv;j++) s->fd_v[j] = 0.; 
  for(i=0;i<nsx;i++){   
	s->fd_x[i] = 0.; 
	for(j=0;j<nsv;j++){
      s->fd_x[i] += s->bfd[i][j];                             /* row fd */
	  s->fd_v[j] += s->bfd[i][j];                             /* col fd */
	}
	s->density_x[i] = s->fd_x[i]/s->num_persons;         /* row density */
  }
  for(j=0;j<nsv;j++) 
	s->density_v[j] = s->fd_v[j]/s->num_persons;         /* col density */

  cum_rel_freqs(0, nsx-1, 1, s->density_x, s->crfd_x);      /* row crfd */
  for (i=0;i<nsx;i++)
    s->prd_x[i] = perc_rank(0, nsx-1, 1, s->crfd_x, static_cast<double>(i));      /* row prd */

  cum_rel_freqs(0, nsv-1, 1, s->density_v, s->crfd_v);      /* col crfd */
  for (j=0;j<nsv;j++)
    s->prd_v[j] = perc_rank(0, nsv-1, 1, s->crfd_v, static_cast<double>(j));      /* col prd */

  /* following code gets cum rel fd as a row-major vector from bfd[][];
     this is needed for bootstrap --- see Parametric_boot_biv() */    

  s->crfd_vector_bfd = dvector(0,ns-1);
  ptr = s->bfd[0];
  for(j=0;j<ns;j++) s->crfd_vector_bfd[j] = *ptr++;     /* assign freqs */   

  for(j=1;j<ns;j++) 
	s->crfd_vector_bfd[j] += s->crfd_vector_bfd[j-1];  /* get cum freqs */
  for(j=0;j<ns;j++) 
	s->crfd_vector_bfd[j] /= s->num_persons;       /* get cum rel freqs */

  /* following code gets rel fd version of bfd[][]; i.e., brfd[][] is
     the smoothed rel freq biv dist for x by v; needed for FEorMFE_EE() */ 

  s->brfd = dmatrix(0,nsx-1,0,nsv-1);
  for(i=0;i<nsx;i++)
    for(j=0;j<nsv;j++) 
	  s->brfd[i][j] = s->bfd[i][j]/s->num_persons;

  return; 
} 

/************************************************************************/

void Smooth_BLL(int n, int nsu, double minu, double incu,
                int nsv, double minv, double incv, double *nct, 
                int anchor, int cu, int cv, int cuv, int cpm[][2], 
				int scale, int Btype, int ctype, double crit,
                FILE *fp, struct BLL_SMOOTH *s)
/*
  Performs bivariate log-linear smoothing in terms of the
  multinomial model as described by Holland and Thayer (1987),
  abbreviated here as H&T.  Assumes v NOT included in u.

  Input
  
    n = number of persons
    nsu = number of score categories for u
    minu = minimum raw score for u
    incu = increment in raw scores for u
    nsv = number of score categories for v
    minv = minimum raw score for v
    incv = increment in raw scores for v
    nct[] = bivariate frequency distribution 
            in row major vector form (see get_nct_bfd())
    cu = number of degrees for smoothing for u
    cv = number of degrees for smoothing for v
    cuv = number of degrees for smoothing for cross-product mts
    cpm[][2] = cross-product moments
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    *fp = pointer to output file
    
  Output: populates struct BLL_SMOOTH s 

  NOTE: nct[] already stored in s by prior call to get_nct_bfd()

  Function calls other than C or NR utilities:
    design_matrix()
    iteration()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,j,
	  max_nit = 40;           /* maximum number of iterations */ 

  int ns = nsu*nsv,    /* # score cats = # rows in design mat */
      nc = cu + cv + cuv;  /* number of cols in design matrix */

  s->num_persons = n;
  s->nsu = nsu;
  s->minu = minu;
  s->incu = incu;
  s->nsv = nsv;
  s->minv = minv;
  s->incv = incv;

  s->anchor = anchor;
  s->cu = cu;
  s->cv = cv;
  s->cuv = cuv;
  s->cpm = imatrix(0,cuv-1,0,1);
  for(i=0;i<cuv;i++) for(j=0;j<2;j++) s->cpm[i][j] = cpm[i][j];

  s->scale = scale;
  s->max_nit = max_nit;
  s->ctype = ctype;
  s->Btype = Btype;
  s->crit = crit;

  s->ns = ns;
  s->nc = nc;

  s->B_raw = dmatrix(0,ns-1,0,nc-1);
  s->B = dmatrix(0,ns-1,0,nc-1);
  s->mct = dvector(0,ns-1);
  s->Beta = dvector(0,nc-1);
  s->n_mts = dvector(0,nc-1);
  s->m_mts = dvector(0,nc-1);
  s->n_mts_raw = dvector(0,nc-1);
  s->m_mts_raw = dvector(0,nc-1);
                      
  design_matrix(nsu,minu,incu, nsv,minv,incv, 
                cu,cv,cuv,cpm, scale,s->B_raw,s->B); 

  /* iteration-step results not printed if first parameter is NULL;
     results are printed if first parameter is fp!=NULL. Note that fp
	 can be set to NULL in Wrapper_Smooth_BLL() */

  s->nit = iteration(fp, s->B, s->B_raw, nct, static_cast<double>(n), 
                     nullptr, ns, cu, cv, cuv, cpm,
                     max_nit, ctype, Btype, crit,
                     s->Beta, s->mct, s->n_mts, s->m_mts, 
                     s->n_mts_raw, s->m_mts_raw, 
                     &(s->lrchisq),&(s->nzero), &(s->ap));     

  return;
}

/********************************************************************************/

void Print_CL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print CL results (CG; log-linear smoothed equipercentile equating)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: 
    score() 
                                            
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
{
  int i,j,
      jlast = (inall->method=='A') ? 75 : 63;

  char x = inall->xv->id1,                            /* x = u[] + v for new form */
       v = inall->xv->id2,                           /* common items for new form */
       y = inall->yv->id1,                            /* y = w[] + z for old form */
       z = inall->yv->id2,                           /* common items for old form */
       u[3],                  /* string designating non-common items for new form */
       w[3];                  /* string designating non-common items for old form */

  if(inall->anchor){u[0] = x; u[1] = '\''; u[2] = '\0';}  /* internal anchor "x'" */
  else {u[0] = ' '; u[1] = x; u[2] = '\0';}               /* external anchor " x" */

  if(inall->anchor){w[0] = y; w[1] = '\''; w[2] = '\0';}  /* internal anchor "y'" */
  else {w[0] = ' '; w[1] = y; w[2] = '\0';}               /* external anchor " y" */
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  std::fprintf(fp,"Equipercentile Equating with CINEG Design"
             " (%s Anchor; w1 = %7.5f)\n"
			 "and Polynomial Log-Linear Smoothing\n\n",
             (inall->anchor != 0) ? "Internal" : "External", inall->w1); 
  std::fprintf(fp,"Input file for %c%c and pop 1: %s\n",x,v,inall->xfname);
  std::fprintf(fp,"Input file for %c%c and pop 2: %s\n\n",y,z,inall->yfname);

  if(inall->anchor)
	  std::fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score; and\n"
				 "(2) variables for old form and pop 2 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score.\n\n",
			     x,u,v,x,u,v,x,u,v,  y,w,z,y,w,z,y,w,z);
  else
	  std::fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
				 "    %c and %c, where %c does not include %c;\n"
				 "(2) variables for old form and pop 2 are identified as\n"
				 "    %c and %c, where %c does not include %c.\n\n",
			     x,v,x,v, y,z,y,z);

  /* information about x and v in pop 1 */

  std::fprintf(fp,"Number of persons for %c = %7d\n\n",x,inall->bllxv->num_persons);
  std::fprintf(fp,"Number of score categories for %c = %d\n",x,inall->bllxv->nsx);
  std::fprintf(fp,"Number of score categories for %c = %d\n",v,inall->bllxv->nsv);

  std::fprintf(fp,"Polynomial degree for %s = %d\n",u,inall->bllxv->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",v,inall->bllxv->cv);
  std::fprintf(fp,"Number of cross-products %s%c = %d\n",u,v,inall->bllxv->cuv);
  if(inall->bllxv->cuv>0) std::fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllxv->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               u,inall->bllxv->cpm[i][0],v,inall->bllxv->cpm[i][1],
                               (inall->bllxv->cuv==i+1) ? ' ' : ',');

 /* information about y and v in pop 2 */

  std::fprintf(fp,"\n\nNumber of persons for %c = %7d\n\n",y,inall->bllyv->num_persons);
  std::fprintf(fp,"Number of score categories for %c = %d\n",y,inall->bllyv->nsx);
  std::fprintf(fp,"Number of score categories for %c = %d\n",z,inall->bllyv->nsv);

  std::fprintf(fp,"Polynomial degree for %s = %d\n",w,inall->bllyv->cu);
  std::fprintf(fp,"Polynomial degree for %c = %d\n",z,inall->bllyv->cv);
  std::fprintf(fp,"Number of cross-products %s%c = %d\n",w,z,inall->bllyv->cuv);
  if(inall->bllyv->cuv>0) std::fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllyv->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               w,inall->bllyv->cpm[i][0],z,inall->bllyv->cpm[i][1],
                               (inall->bllyv->cuv==i+1) ? ' ' : ',');
  
  /* print equated raw scores for methods */

  if(inall->method == 'E' || inall->method == 'G' || inall->method == 'A')
    std::fprintf(fp,"\n\nBraun-Holland (BH) Under Frequency Estimation (FE):\n"
               "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
               r->b[0], r->a[0]);
  if(inall->method == 'F' || inall->method == 'G' || inall->method == 'A')
    std::fprintf(fp,"Braun-Holland (BH) Under Modified Frequency Estimation (MFE):\n"
               "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
               r->b[1], r->a[1]);

  if(inall->method == 'E' || inall->method == 'F' || 
	 inall->method == 'G' || inall->method == 'A')
     for(i=1;i<=jlast;i++) std::fprintf(fp,"-");

  /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (%c)  ",x);
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"  %s",inall->names[j]);
  std::fprintf(fp,"\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    std::fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  std::fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][0]);
  std::fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][1]);
  std::fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][2]);
  std::fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][3]);
  
  std::fprintf(fp,"\n\n");
  for(j=1;j<=jlast;j++) std::fprintf(fp,"*");

}

/****************************************************************************/

void Wrapper_Smooth_ULL(struct USTATS *x, int c,
						int scale, int Btype, int ctype, double crit,
						FILE *fp, struct ULL_SMOOTH *s)
/*
  Wrapper to do univariate log-linear smoothing.

  Input
    x   =  UTATS structure
    c   = number of degrees for polynomial smoothing
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    fp  = pointer to output file (if NULL then no output printed;
	      in particular, results for each iteration step are not printed)

  Output
    populates s, which is a ULL_SMOOTH structure

  Function calls other than C or NR utilities:  
    Smooth_ULL()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  Smooth_ULL(x->n, x->ns, x->min, x->inc, x->dbl_fd, c, 
	         scale, Btype, ctype, crit, fp, s); 

  return; 
} 

/**************************************************************/

void Smooth_ULL(int n, int ns, double min, double inc,
                double *fd, int c, 
				int scale, int Btype, int ctype, double crit,
				FILE *fp, struct ULL_SMOOTH *s)
/*
  Performs univariate log-linear smoothing in terms of the 
  multinomial model as described by Holland and Thayer (1987),
  abbreviated here as H&T.  

  Input
  
    n = number of persons
    ns = number of score categories
    min = minimum raw score
    inc = increment in raw scores
    fd[] = frequency distribution
    c = number of degrees for polynomial smoothing
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    *fp = pointer to output file
    
  Output: populates struct ULL_SMOOTH s 

  Function calls other than C or NR utilities:
    design_matrix()
    iteration()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,
      max_nit = 40;           /* maximum number of iterations */ 

  s->num_persons = n;
  s->ns = ns;
  s->min = min;
  s->inc = inc;
  s->c = c;
  s->nct = fd;
  s->scale = scale;
  s->max_nit = max_nit;
  s->ctype = ctype;
  s->Btype = Btype;
  s->crit = crit;

  s->B_raw = dmatrix(0,ns-1,0,c-1);
  s->B = dmatrix(0,ns-1,0,c-1);
  s->mct = dvector(0,ns-1);
  s->Beta = dvector(0,c-1);
  s->n_mts = dvector(0,c-1);
  s->m_mts = dvector(0,c-1);
  s->n_mts_raw = dvector(0,c-1);
  s->m_mts_raw = dvector(0,c-1);

  s->density = dvector(0,ns-1);
  s->crfd = dvector(0,ns-1);
  s->prd = dvector(0,ns-1);
                      
  design_matrix(ns,min,inc, 0,0,0, c,0,0,nullptr, scale,s->B_raw,s->B); 

  /* iteration-step results not printed if first parameter is NULL;
     results are printed if first parameter is fp!=NULL. Note that fp
	 can be set to NULL in Wrapper_Smooth_ULL() */ 

  s->nit = iteration(fp, s->B, s->B_raw, fd, static_cast<double>(n), 
                     nullptr, ns, c, 0, 0, nullptr,
                     max_nit, ctype, Btype, crit,
                     s->Beta, s->mct, s->n_mts, s->m_mts, 
                     s->n_mts_raw, s->m_mts_raw, 
                     &(s->lrchisq),&(s->nzero),&(s->ap));

  for (i=0;i<ns;i++) s->density[i] = s->mct[i]/s->num_persons;
  cum_rel_freqs(0, ns-1, 1, s->density, s->crfd); 
  for (i=0;i<ns;i++)
    s->prd[i] = perc_rank(0, ns-1, 1, s->crfd, static_cast<double>(i));

  return;
}

/*******************************************************************/

void Print_ULL(FILE *fp, char tt[], struct USTATS *x,
               struct ULL_SMOOTH *s, int print_dm, int print_mts)
/*
  print results in ULL_SMOOTH s
  
  Input
    fp   = file pointer for output
    tt[] = user supplied text identifier
    x    = struct USTATS 
    s    = struct ULL_SMOOTH
    print_dm =  print design matrices
	            0 --> no printing of matrices 
	            1 --> print B (and coefficients of B)
                2 --> print B and B_raw (and coefficients of B)
				NOTE: if scale = 0, then only B_raw is printed
   print_mts = print moments
	            0 --> no print
				1 --> print moments for B
				2 --> print moments for B and B-raw
				NOTE: if scale = 0, then only moments of B_raw are printed
                   					
  Function calls other than C or NR utilities: 
    score()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j;

  std::fprintf(fp,"");

  std::fprintf(fp,"\n\n%s\n\n",tt);
  std::fprintf(fp,"Input filename:  %s\n\n",x->fname);

  std::fprintf(fp,"Univariate Log-Linear Smoothing using the Multinomial Model\n\n");

  std::fprintf(fp,"Number of score categories = %d\n",s->ns);
  std::fprintf(fp,"Number of persons (i.e., total of frequencies) = %7d\n\n",
	         s->num_persons);

  std::fprintf(fp,"Polynomial degree = %d\n\n",s->c);

  if(s->scale==0) std::fprintf(fp,"Design matrix B is powers of raw scores\n\n");
  else std::fprintf(fp,"Design matrix B is scaled powers of raw scores\n"
                  "such that for any column of B, sum(elements) = 0 and\n"
                  "sum(elements^2) = 1\n\n");
  
  std::fprintf(fp,"Convergence criterion:  ");
  if(s->Btype==0 && s->ctype==0)
    std::fprintf(fp,"|(n_mts[i] - m_mts[i]| <= %15.10f",s->crit); 
  else if(s->Btype==0 && s->ctype==1)
    std::fprintf(fp,"|n_mts[i] - m_mts[i])/n_mts[i]| <= %15.10f",s->crit);
  else if(s->Btype==1 && s->ctype==0)
    std::fprintf(fp,"|(n_mts_raw[i] - m_mts_raw[i]| <= %15.10f",s->crit); 
  else 
    std::fprintf(fp,"|n_mts_raw[i] - m_mts_raw[i])/n_mts_raw[i]|"
                 " <= %15.10f",s->crit);
  std::fprintf(fp,"  for i = 1 to %d\n\n",s->c);

  std::fprintf(fp,"Number of iterations to convergence = %d\n\n",s->nit);

  std::fprintf(fp,"Likelihood-ratio chi-square = %12.5f "
             "with %d degrees of freedom\n\n",s->lrchisq,s->ns-s->c-1);

  /* NOTE: if scale=0, the B_raw matrix is the design matric B.  In 
     this case only B_raw is printed along with its coefficients
	 and mts_raw[] */

  if(s->scale==0) std::fprintf(fp,"Since scale=0, B_raw=B is the design matrix\n\n");

  /* print coefficients for design matrix */

  std::fprintf(fp,"COEFFICIENTS:\n\n");
    for(i=1;i<=s->c;i++)
      std::fprintf(fp,"   Beta[%2d] = %12.5f\n",i,s->Beta[i-1]); 

  /* print B matrix (scaled) if print_dm is positive and scale==1 */

  if(print_dm && s->scale==1){ 

    std::fprintf(fp,"\n\nB DESIGN MATRIX");
    std::fprintf(fp,"\n\nCategory   Score");
    for(j=0;j<s->c;j++) std::fprintf(fp,"    Deg %2d",j+1); 
    std::fprintf(fp,"\n"); 
    for(i=0;i<s->ns;i++){
      std::fprintf(fp,"\n%8d%8.3f",i,score(i,s->min,s->inc));
      for(j=0;j<s->c;j++)
        std::fprintf(fp,"%10.5f",s->B[i][j]);
    }
  }

  /* print B_raw matrix if print_dm==2 and scale==1,
     or print_dm==1 and scale==0 */

  if( (print_dm==2 && s->scale==1)  || (print_dm && s->scale==0) ){

    std::fprintf(fp,"\n\nB_raw DESIGN MATRIX");
    std::fprintf(fp,"\n\nCategory   Score");
    for(j=0;j<s->c;j++) std::fprintf(fp,"           Degree %2d",j+1);
    std::fprintf(fp,"\n"); 
    for(i=0;i<s->ns;i++){
      std::fprintf(fp,"\n%8d%8.3f",i,score(i,s->min,s->inc));
      for(j=0;j<s->c;j++)
        std::fprintf(fp,"%20.5f",s->B_raw[i][j]);
    }
  }

  /* print actual and fitted frequencies */

  std::fprintf(fp,"\n\n");
  for(i=1;i<=93;i++) std::fprintf(fp,"-");

  std::fprintf(fp,"\n\n                      Actual      Fitted\n");
  std::fprintf(fp,"                      freqs:      freqs:\n");
  std::fprintf(fp,"Category   Score       nct[]       mct[]          prop  fitted-prop"
	         "  fitted-crfd   fitted-prd\n\n");
  for(i=0;i<s->ns;i++)
    std::fprintf(fp,"%8d%8.3f%12.0f%12.5f       %7.5f %12.5f %12.5f %12.5f\n",
	        i,score(i,s->min,s->inc),s->nct[i],s->mct[i],
			s->nct[i]/s->num_persons,s->density[i],
			s->crfd[i],s->prd[i]);

  /* print moments of B matrix (scaled) if print_mts is positive and scale==1 */

  if(print_mts && s->scale==1){
    std::fprintf(fp,"\n");
    for(i=1;i<=s->c;i++)
      std::fprintf(fp,"         mts[%2d]%12.5f%12.5f\n",
                i,s->n_mts[i-1], s->m_mts[i-1]);  
  }

  /* print B_raw moments if print_mts==2 and scale==1,
     or print_mts==1 and scale==0 */

  if( (print_mts==2 && s->scale==1)  || (print_mts && s->scale==0) ){
    std::fprintf(fp,"\n");
    for(i=1;i<=s->c;i++)
      std::fprintf(fp,"     mts_raw[%2d]%12.5f%12.5f\n",
                 i,s->n_mts_raw[i-1], s->m_mts_raw[i-1]);  
  }

  if(print_mts) std::fprintf(fp,"\nNOTES:\n\n");

  if(print_mts && s->scale==1){	
    std::fprintf(fp,  
      " mts[] are moments associated with the design matrix B that is\n"
      " used to obtain the log-linear solution. Columns of B are either:\n" 
      " (a) powers of raw scores (using min and inc); this is\n"  
      "     accomplished by setting scale=0 in Smooth_BLL(); or\n" 
      " (b) scaled powers of raw scores such that for each column of\n"  
      "     B, sum(elements) = 0 and sum(elements^2) = 1; this is\n"
      "     accomplished by setting scale=1 in Smooth_BLL().\n" 
      " (Experience suggests that scale=1 leads to more numerical\n"    /* edit change 3-19-09 */
      " stability.) Letting i be a score category, the jth moment is\n" /* edit change 3-19-09 */
      " mts[j] = sum_i(B[i][j-1]*rf[i]) where rf[] is either actual\n" 
      " relative frequencies or fitted relative frequencies.\n\n");
  }

   if( (print_mts==2 && s->scale==1)  || (print_mts && s->scale==0) ){
    std::fprintf(fp,
      " mts_raw[] are based on (a) above (i.e, no scaling). Importantly,\n" 
      " however, the moments greater than or equal to 3 are central moments\n" 
      " associated with matrix B_raw in the following sense. Letting i be a\n"
      " score category, rf[i] be a relative frequency, and mean = sum_i(B_raw[i][0]*rf[i]),\n"
	  " then:\n"
      " mts_raw[1] = mean,  mts_raw[2] = sd, and for j = 3,...,%d,\n" 
      " mts_raw[j] =  sum_i(B_raw[i][0] - mean)^(j)*rf[i]/sd^(j).\n\n",s->c);
  }
  
  for(i=1;i<=70;i++)  std::fprintf(fp,"*");

  return;
} 

/*******************************************************************/ 

void Wrapper_Smooth_BLL(struct BSTATS *xv, int anchor,
                        int cu, int cv, int cuv, int cpm[][2],
						int scale, int Btype, int ctype, double crit,
                        FILE *fp, struct BLL_SMOOTH *s)
/*
  Wrapper to do univariate log-linear smoothing.

  NOTE: if this function is used with the SG design,
        anchor should be set to 0 (external)

  Input
    xv      =  BTATS structure
    anchor  = 0 (external); = 1 (internal) 
    cu      = number of degrees of smoothing for u
    cv      = number of degrees of smoothing for v
    cuv     = number of cross-product moments
    cpm[cuv-1][2] = zero-offset matrix designating cross-
                    product moments.  Example: let cuv = 3,
                    and the desired cross-product moments be 
                    (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                    Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                    cpm[2] = {2,1}. 
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    fp  = pointer to output file (if NULL then no output printed;
	      in particular, results for each iteration step are not printed) 

  Output
    gets s->nct and 
    calls function that populates other elements of s

  NOTE: When anchor is internal, smoothing is done on matrix of
       (non-common-item scores) by (common-item scores).  For an
	   internal anchor, smoothing could be done on matrix of
	   (total scores) by (common-item scores), and convergence
	   should not be a problem, but structural zeros likely would be
	   given fitted positive values.  (For an internal anchor,
	   structural zeros occur when a total score is impossible
	   given a particular common item score.)

  Function calls other than C or NR utilities:  
    get_nct_bfd()
    Smooth_BLL()
    get_bfd_mct()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j,
      nsx = xv->ns1, 
      nsv = xv->ns2,
	  ns = nsx*nsv,
      nsu = (anchor) ? xv->ns1 - xv->ns2 + 1: xv->ns1;
  double *ptr,                                /* pointer variable */
	     minu = (anchor) ?  xv->min1 - xv->min2: xv->min1,
         incu = xv->inc1;    /* to be sensible, it must be true that
                                       xv->inc1 = xv->inc2 = incu */

  /* nsu is number of categories for non-common items */
 
  s->nct = dvector(0,nsu*nsv-1);
  s->bfd = dmatrix(0,nsx-1,0,nsv-1);  

  /* For an external anchor, directly convert xv->bfd[][] 
     to a row-major vector nct[]. For an internal anchor, 
     conceptually we first collapse xv->bfd[][] to 
     uv->bfd[][], and then convert to nct[].  
     In the end this means that nct[] is a row-major vector 
     with nsu*nsv elements */

  get_nct_bfd(anchor, xv->ns1, xv->ns2, xv->dbl_bfd, s->nct); 

  /* call bivariate log-linear smoothing function, which
     assumes that row categores are NOT included in
     column categories */ 

  Smooth_BLL(xv->n, nsu, minu, incu, 
             xv->ns2, xv->min2, xv->inc2, s->nct, 
             anchor, cu, cv, cuv, cpm, 
			 scale, Btype, ctype, crit, fp, s); 

  /* For an external anchor, directly convert the row major 
     vector s->mct[nsu*nsv] to s->bfd[nsu][nsv]. For an 
     internal anchor, convert s->mct[nsu*nsv] to
     s->bfd[nsu+nsv-1][nsv] repositioning elements and
	 adding structural zeros. */

  get_bfd_mct(anchor, xv->ns1, xv->ns2, s->mct, s->bfd);

  /* get row and col marginal densities, crfd's and prd's for bfd[][] */

  s->nsx = nsx;
  s->minx = xv->min1;
  s->incx = xv->inc1;

  s->fd_x = dvector(0,nsx-1);    /* row marginal frequencies for bfd[][] */
  s->density_x = dvector(0,nsx-1);   /* row marginal density for bfd[][] */
  s->crfd_x = dvector(0,nsx-1);         /* row marginal crfd for bfd[][] */
  s->prd_x = dvector(0,nsx-1);            /* row marginal PR for bfd[][] */

  s->fd_v = dvector(0,nsv-1);    /* col marginal frequencies for bfd[][] */
  s->density_v = dvector(0,nsv-1);   /* col marginal density for bfd[][] */
  s->crfd_v = dvector(0,nsv-1);         /* col marginal crfd for bfd[][] */
  s->prd_v = dvector(0,nsv-1);            /* col marginal PR for bfd[][] */
                      
  for(j=0;j<nsv;j++) s->fd_v[j] = 0.; 
  for(i=0;i<nsx;i++){   
	s->fd_x[i] = 0.; 
	for(j=0;j<nsv;j++){
      s->fd_x[i] += s->bfd[i][j];                             /* row fd */
	  s->fd_v[j] += s->bfd[i][j];                             /* col fd */
	}
	s->density_x[i] = s->fd_x[i]/s->num_persons;         /* row density */
  }
  for(j=0;j<nsv;j++) 
	s->density_v[j] = s->fd_v[j]/s->num_persons;         /* col density */

  cum_rel_freqs(0, nsx-1, 1, s->density_x, s->crfd_x);      /* row crfd */
  for (i=0;i<nsx;i++)
    s->prd_x[i] = perc_rank(0, nsx-1, 1, s->crfd_x, static_cast<double>(i));      /* row prd */

  cum_rel_freqs(0, nsv-1, 1, s->density_v, s->crfd_v);      /* col crfd */
  for (j=0;j<nsv;j++)
    s->prd_v[j] = perc_rank(0, nsv-1, 1, s->crfd_v, static_cast<double>(j));      /* col prd */

  /* following code gets cum rel fd as a row-major vector from bfd[][];
     this is needed for bootstrap --- see Parametric_boot_biv() */    

  s->crfd_vector_bfd = dvector(0,ns-1);
  ptr = s->bfd[0];
  for(j=0;j<ns;j++) s->crfd_vector_bfd[j] = *ptr++;     /* assign freqs */   

  for(j=1;j<ns;j++) 
	s->crfd_vector_bfd[j] += s->crfd_vector_bfd[j-1];  /* get cum freqs */
  for(j=0;j<ns;j++) 
	s->crfd_vector_bfd[j] /= s->num_persons;       /* get cum rel freqs */

  /* following code gets rel fd version of bfd[][]; i.e., brfd[][] is
     the smoothed rel freq biv dist for x by v; needed for FEorMFE_EE() */ 

  s->brfd = dmatrix(0,nsx-1,0,nsv-1);
  for(i=0;i<nsx;i++)
    for(j=0;j<nsv;j++) 
	  s->brfd[i][j] = s->bfd[i][j]/s->num_persons;

  return; 
} 

/************************************************************************/

void Smooth_BLL(int n, int nsu, double minu, double incu,
                int nsv, double minv, double incv, double *nct, 
                int anchor, int cu, int cv, int cuv, int cpm[][2], 
				int scale, int Btype, int ctype, double crit,
                FILE *fp, struct BLL_SMOOTH *s)
/*
  Performs bivariate log-linear smoothing in terms of the
  multinomial model as described by Holland and Thayer (1987),
  abbreviated here as H&T.  Assumes v NOT included in u.

  Input
  
    n = number of persons
    nsu = number of score categories for u
    minu = minimum raw score for u
    incu = increment in raw scores for u
    nsv = number of score categories for v
    minv = minimum raw score for v
    incv = increment in raw scores for v
    nct[] = bivariate frequency distribution 
            in row major vector form (see get_nct_bfd())
    cu = number of degrees for smoothing for u
    cv = number of degrees for smoothing for v
    cuv = number of degrees for smoothing for cross-product mts
    cpm[][2] = cross-product moments
    scale = type of scaling:
	        0 --> no scaling; 
            1 --> scale such that each column of B has
                  sum (elements) = 0 and sum (elements^2) = 1
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
	ctype = comparison type for criterion:
	        0 --> means use absolute criterion; 
			1 --> means use relative criterion
	crit = convergence criterion value
    *fp = pointer to output file
    
  Output: populates struct BLL_SMOOTH s 

  NOTE: nct[] already stored in s by prior call to get_nct_bfd()

  Function calls other than C or NR utilities:
    design_matrix()
    iteration()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,j,
	  max_nit = 40;           /* maximum number of iterations */ 

  int ns = nsu*nsv,    /* # score cats = # rows in design mat */
      nc = cu + cv + cuv;  /* number of cols in design matrix */

  s->num_persons = n;
  s->nsu = nsu;
  s->minu = minu;
  s->incu = incu;
  s->nsv = nsv;
  s->minv = minv;
  s->incv = incv;

  s->anchor = anchor;
  s->cu = cu;
  s->cv = cv;
  s->cuv = cuv;
  s->cpm = imatrix(0,cuv-1,0,1);
  for(i=0;i<cuv;i++) for(j=0;j<2;j++) s->cpm[i][j] = cpm[i][j];

  s->scale = scale;
  s->max_nit = max_nit;
  s->ctype = ctype;
  s->Btype = Btype;
  s->crit = crit;

  s->ns = ns;
  s->nc = nc;

  s->B_raw = dmatrix(0,ns-1,0,nc-1);
  s->B = dmatrix(0,ns-1,0,nc-1);
  s->mct = dvector(0,ns-1);
  s->Beta = dvector(0,nc-1);
  s->n_mts = dvector(0,nc-1);
  s->m_mts = dvector(0,nc-1);
  s->n_mts_raw = dvector(0,nc-1);
  s->m_mts_raw = dvector(0,nc-1);
                      
  design_matrix(nsu,minu,incu, nsv,minv,incv, 
                cu,cv,cuv,cpm, scale,s->B_raw,s->B); 

  /* iteration-step results not printed if first parameter is NULL;
     results are printed if first parameter is fp!=NULL. Note that fp
	 can be set to NULL in Wrapper_Smooth_BLL() */

  s->nit = iteration(fp, s->B, s->B_raw, nct, static_cast<double>(n), 
                     nullptr, ns, cu, cv, cuv, cpm,
                     max_nit, ctype, Btype, crit,
                     s->Beta, s->mct, s->n_mts, s->m_mts, 
                     s->n_mts_raw, s->m_mts_raw, 
                     &(s->lrchisq),&(s->nzero), &(s->ap));     

  return;
}

/********************************************************************************/

void Print_BLL(FILE *fp, char tt[], struct BSTATS *xv, struct BLL_SMOOTH *s, 
			   int print_dm, int print_mts, int print_freq, int print_bfd)
/*
  print results in BLL_SMOOTH s
  
  Input
    fp   = file pointer for output
    tt[] = user supplied text identifier
    xv   = struct BSTATS 
    s    = struct BLL_SMOOTH
    print_dm = print design matrices:
                 0 --> no printing  
                 1 --> print B (and coefficients for B);
                 2 --> print B and B_raw (and coefficients for B)
    print_mts = print moments
	             0 --> no print
				 1 --> print moments for B
				 2 --> print moments for B and B-raw
    print_freq = format for printing actual and fitted frequencies
	             for (non-common) by (common)  
	             0 --> don't print
				 1 --> column vector format 
                 2 --> matrix format
    print_bfd = print bfd[][] (nsx rows by nsv cols where nsx = nsu 
	             if anchor is external)
				 0 --> don't print bfd and marginals
				 1 --> print bfd and marginals
				 2 --> print marginals only

    NOTE: For an external anchor print_freq==1 and print_bfd==1 
	      cause the same matrix to be printed twice

  Function calls other than C or NR utilities: 
    score()

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j,cells;
  char x = xv->id1;
  char v = xv->id2;
  char u[3];                                 /* string designating non-common items */
  if(s->anchor){u[0] = x; u[1] = '\''; u[2] = '\0';}  /* int anch: x' for non-common*/
  else {u[0] = ' '; u[1] = x; u[2] = '\0';}               /* external anchor x = x' */
 
  std::fprintf(fp,"");

  std::fprintf(fp,"\n\n%s\n\n",tt);
  std::fprintf(fp,"Input filename:  %s\n\n",xv->fname);

  std::fprintf(fp,"Bivariate Log-Linear Smoothing using the Multinomial Model\n\n");

  if(s->anchor)
	std::fprintf(fp,"NOTE: In this output, variables are identified as\n"
	           "%c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	           "total score, %s is the non-common-item score, and\n"
			   "%c is the common-item score.\n\n",
			   x,u,v,x,u,v,x,u,v);
  else
    std::fprintf(fp,"NOTE: In this output, variables are identified as %c and %c.\n"
	           "Since the anchor is external (explcitly if the design is CG;\n"
			   "implicitly if the design is SG), %c does not include %c\n\n",
			   x,v,x,v);
 
  std::fprintf(fp,"Number of score categories for %s  = %d\n",u,s->nsu);
  std::fprintf(fp,"Number of score categories for %c  = %d\n",v,s->nsv);
  std::fprintf(fp,"Total number of score categories = %d\n\n",s->ns);

  std::fprintf(fp,"Number of persons (i.e., total of frequencies) = %7d\n\n",
	         s->num_persons);

  std::fprintf(fp,"Polynomial degree for %s  = %d\n",u,s->cu);
  std::fprintf(fp,"Polynomial degree for %c  = %d\n",v,s->cv);
  std::fprintf(fp,"Number of cross-products (%s,%c) = %d\n",u,v,s->cuv);
  if(s->cuv>0){
    std::fprintf(fp,"Cross-Product moments: ");
    for(i=0;i<s->cuv;i++) std::fprintf(fp,"  (%s%d,%c%d)%c",
                               u,s->cpm[i][0],v,s->cpm[i][1],
                               (s->cuv==i+1) ? ' ' : ',');
  }
  std::fprintf(fp,"\nNumber of columns in design matrix = %d\n\n",s->nc);

  if(s->scale==0) std::fprintf(fp,"Design matrix B is powers of raw scores\n\n");
  else std::fprintf(fp,"Design matrix B is scaled powers of raw scores\n"
                  "such that for any column of B, sum(elements) = 0 and\n"
                  "sum(elements^2) = 1\n\n");
  
  std::fprintf(fp,"Convergence criterion:  ");
  if(s->Btype==0 && s->ctype==0)
    std::fprintf(fp,"|(n_mts[i] - m_mts[i]| <= %15.10f",s->crit); 
  else if(s->Btype==0 && s->ctype==1)
    std::fprintf(fp,"|n_mts[i] - m_mts[i])/n_mts[i]| <= %15.10f",s->crit);
  else if(s->Btype==1 && s->ctype==0)
    std::fprintf(fp,"|(n_mts_raw[i] - m_mts_raw[i]| <= %15.10f",s->crit); 
  else 
    std::fprintf(fp,"|n_mts_raw[i] - m_mts_raw[i])/n_mts_raw[i]|"
                 " <= %15.10f",s->crit);
  std::fprintf(fp,"  for i = 1 to %d\n\n",s->nc);

  std::fprintf(fp,"Number of iterations to convergence = %d\n\n",s->nit);

  std::fprintf(fp,"Likelihood-ratio chi-square = %12.5f "
             "with %d degrees of freedom\n\n",s->lrchisq,
              s->ns  - s->nc - 1);

  if(s->scale==0) std::fprintf(fp,"Since scale=0, B_raw=B is the design matrix\n\n");

  /* print coefficients */

  
  std::fprintf(fp,"COEFFICIENTS:\n\n");
    for(i=1;i<=s->nc;i++)
      std::fprintf(fp,"   Beta[%2d] = %12.5f\n",i,s->Beta[i-1]);
  
   /* print B matrix (scaled) if print_dm is positive and scale==1 */

  if(print_dm && s->scale==1){

    std::fprintf(fp,"\n\nB DESIGN MATRIX");
    std::fprintf(fp,"\n\nCategory");
    for(j=0;j<s->nc;j++) std::fprintf(fp,"    Deg %2d",j+1); 
    std::fprintf(fp,"\n"); 
    for(i=0;i<s->ns;i++){
      std::fprintf(fp,"\n%8d",i);
      for(j=0;j<s->nc;j++)
        std::fprintf(fp,"%10.5f",s->B[i][j]);
    }
  }

  /* print B_raw matrix if print_dm==2 and scale==1,
     or print_dm==1 and scale==0 */

  if( (print_dm==2 && s->scale==1)  || (print_dm && s->scale==0) ){

    std::fprintf(fp,"\n\nB_raw DESIGN MATRIX");
    std::fprintf(fp,"\n\nCategory");
    for(j=0;j<s->nc;j++) std::fprintf(fp,"           Degree %2d",j+1);
    std::fprintf(fp,"\n"); 
    for(i=0;i<s->ns;i++){
      std::fprintf(fp,"\n%8d",i);
      for(j=0;j<s->nc;j++)
        std::fprintf(fp,"%20.5f",s->B_raw[i][j]);
    }
  }

  /* print actual and fitted frequencies */

  if(print_freq==1){                                /* vector format */
    std::fprintf(fp,"\n\n                      Actual      Fitted\n");
    std::fprintf(fp,"                      freqs:      freqs:\n");
    std::fprintf(fp,"Category               nct[]       mct[]\n\n");
    for(i=0;i<s->ns;i++)
      std::fprintf(fp,"%8d        %12.0f%12.5f\n",i,
              s->nct[i],s->mct[i]);
  }

  if(print_freq==2){                                /* matrix format */
    std::fprintf(fp,"\n\n\nACTUAL FREQUENCIES\n");

    std::fprintf(fp,"\nCategory        ");
    for(j=0;j<s->nsv;j++) std::fprintf(fp,"   %c=%3d",v,j);
    std::fprintf(fp,"\n           Score");
    for(j=0;j<s->nsv;j++) std::fprintf(fp,"%8.3f",score(j,s->minv,s->incv));
    std::fprintf(fp,"\n");

    cells = 0;
    for(i=0;i<s->nsu;i++){
      std::fprintf(fp,"\n  %s=%3d%8.3f",u,i,score(i,s->minu,s->incu));
      for(j=0;j<s->nsv;j++)
        std::fprintf(fp,"%8.0f",s->nct[cells++]);
    }

    std::fprintf(fp,"\n\n\nFITTED FREQUENCIES\n");

    std::fprintf(fp,"\nCategory        ");
    for(j=0;j<s->nsv;j++) std::fprintf(fp,"   %c=%3d",v,j);
    std::fprintf(fp,"\n           Score");
    for(j=0;j<s->nsv;j++) std::fprintf(fp,"%8.3f",score(j,s->minv,s->incv));
    std::fprintf(fp,"\n\n");

    cells = 0;
    for(i=0;i<s->nsu;i++){
      std::fprintf(fp,"\n  %s=%3d%8.3f",u,i,score(i,s->minu,s->incu));
      for(j=0;j<s->nsv;j++)
        std::fprintf(fp,"%8.3f",s->mct[cells++]);
    }
  }

  /* print moments of B matrix (scaled) if print_mts is positive and scale==1 */

  if(print_mts && s->scale==1){

    std::fprintf(fp,"\n\nMOMENTS\n");
    std::fprintf(fp,"                    Based on    Based on\n");
    std::fprintf(fp,"                      Actual      Fitted\n");
    std::fprintf(fp,"                       nct[]       mct[]");

    std::fprintf(fp,"\n\n%s moments based on B:\n\n",u);
    for(i=1;i<=s->cu;i++)
      std::fprintf(fp,"         mts[%2d]%12.5f%12.5f\n",
                i,s->n_mts[i-1], s->m_mts[i-1]);  

    std::fprintf(fp,"\n %c moments based on B:\n\n",v);
    for(i=1;i<=s->cv;i++)
      std::fprintf(fp,"         mts[%2d]%12.5f%12.5f\n",
                 i,s->n_mts[s->cu+i-1], s->m_mts[s->cu+i-1]);  

    std::fprintf(fp,"\n%s%c moments based on B:\n\n",u,v);
    for(i=1;i<=s->cuv;i++)
      std::fprintf(fp,"         mts[%2d]%12.5f%12.5f\n",
                 i,s->n_mts[s->cu + s->cv + i - 1], 
                   s->m_mts[s->cu +s ->cv + i - 1]);
  }

  /* print B_raw mmoments if print_mts==2 and scale==1,
     or print_mts==1 and scale==0 */

  if( (print_mts==2 && s->scale==1)  || (print_mts && s->scale==0) ){

    std::fprintf(fp,"\n\n%s (central) moments based on B_raw:\n\n",u);
    for(i=1;i<=s->cu;i++)
      std::fprintf(fp,"     mts_raw[%2d]%12.5f%12.5f\n",
                i,s->n_mts_raw[i-1], s->m_mts_raw[i-1]);  

    std::fprintf(fp,"\n %c (central) moments based on B_raw:\n\n",v);
    for(i=1;i<=s->cv;i++)
      std::fprintf(fp,"     mts_raw[%2d]%12.5f%12.5f\n",
                 i,s->n_mts_raw[s->cu+i-1], s->m_mts_raw[s->cu+i-1]);  

    std::fprintf(fp,"\n%s%c (central) moments based on B_raw:\n\n",u,v);
    for(i=1;i<=s->cuv;i++)
      std::fprintf(fp,"     mts_raw[%2d]%12.5f%12.5f\n",
                 i,s->n_mts_raw[s->cu + s->cv + i - 1], 
                   s->m_mts_raw[s->cu +s ->cv + i - 1]);
  }

  if(print_mts) std::fprintf(fp,"\nNOTES:\n\n");

  if(print_mts && s->scale==1){	
    std::fprintf(fp,
      " mts[] are moments associated with the design matrix B that is\n"
      " used to obtain the log-linear solution. Columns of B are either:\n" 
      " (a) powers of raw scores (using min and inc); this is\n"  
      "     accomplished by setting scale=0 in Smooth_BLL(); or\n" 
      " (b) scaled powers of raw scores such that for each column of\n"  
      "     B, sum(elements) = 0 and sum(elements^2) = 1; this is\n"
      "     accomplished by setting scale=1 in Smooth_BLL().\n" 
      " (Experience suggests that scale=1 is more likely to lead to\n" 
      " convergence.) Letting i be a score category, the jth moment is\n" 
      " mts[j] = sum_i(B[i][j-1]*rf[i]) where rf[] is either actual\n" 
      " relative frequencies or fitted relative frequencies.\n\n");
  }

   if( (print_mts==2 && s->scale==1)  || (print_mts && s->scale==0) ){
    std::fprintf(fp,
      " mts_raw[] are based on (a) above (i.e, no scaling). Importantly,\n" 
      " however, the moments greater than or equal to 3 are central moments\n" 
      " associated with matrix B_raw in the following sense.\n"
      " (i)   The first %d columns of B_raw are associated with the %d moments\n"
      "       of %s:  mts_raw[1] = mean,  mts_raw[2] = sd, and for j = 3,...,%d,\n" 
      "       mts_raw[j] =  sum_i(B_raw[i][j-1] - mean)^(j)*rf[i]/sd^(j),\n" 
      "       where mean = sum_i(B_raw[i][0]*f[i]).\n"
      " (ii)  The next %d columns of B_raw are associated with the\n" 
      "       (central) moments of of %c.\n",
      s->cu, s->cu, u, s->cu,    s->cv, v);
    if(s->cuv==1)  std::fprintf(fp,
      " (iii) The last columns of B_raw is associated with a (central)\n"
      "       cross-product moment.  (For example, if cpm[0][0] = 1 and\n"
      "       cpm[0][1] = 1, then the cross-product moment is\n"  
      "       the correlation between %s and %c.)\n\n",u,v);
    if(s->cuv>1)  std::fprintf(fp,
	  " (iii) The next %d columns of B_raw are associated with %d (central)\n"
      "       cross-product moments.  (For example, if cpm[0][0] = 1 and\n"
      "       cpm[0][1] = 1, then the 1st cross-product moment is\n"  
      "       the correlation between %s and %c.)\n\n",s->cuv, s->cuv,u,v);
  }

  for(i=1;i<=70;i++)  std::fprintf(fp,"*");

  if(print_bfd==1){                                          /* print bfd[][] */
                                      
	std::fprintf(fp,"\n\nFITTED BIVARIATE FREQUENCY DISTRIBUTION OF %c AND %c: "
		       "bfd[][]\n\n",x,v);
	std::fprintf(fp,"Category        ");
    for(j=0;j<s->nsv;j++) std::fprintf(fp,"   %c=%3d",v,j);
    std::fprintf(fp,"\n           Score");
    for(j=0;j<s->nsv;j++) std::fprintf(fp,"%8.3f",score(j,s->minv,s->incv));
    std::fprintf(fp,"\n\n");

    for(i=0;i<s->nsx;i++){
      std::fprintf(fp,"\n   %c=%3d%8.3f",x,i,score(i,s->minu,s->incu));
      for(j=0;j<s->nsv;j++) std::fprintf(fp,"%8.3f",s->bfd[i][j]);
	}
  }
    
  if(print_bfd){  /* print row and col marginal freq, density, crfd, and prd */

	std::fprintf(fp,"\n\n\n                         "
		       "FITTED MARGINAL DISTRIBUTION OF %c\n",x);  /* print for rows */
    for(i=1;i<=17;i++) std::fprintf(fp," ");
	for(i=1;i<=51;i++) std::fprintf(fp,"-");
    std::fprintf(fp,"\nCategory   Score         Freq     Rel Freq    Cum RFreq"
	         "    Perc Rank\n\n");
    for(i=0;i<s->nsx;i++)
      std::fprintf(fp,"%8d%8.3f %12.5f %12.5f %12.5f %12.5f\n",
	          i,score(i,s->minx,s->incx),s->fd_x[i],s->density_x[i],
	          s->crfd_x[i],s->prd_x[i]);

	std::fprintf(fp,"\n\n                         "
		       "FITTED MARGINAL DISTRIBUTION OF %c\n",v);  /* print for cols */
    for(i=1;i<=17;i++) std::fprintf(fp," ");
    for(i=1;i<=51;i++) std::fprintf(fp,"-");
    std::fprintf(fp,"\nCategory   Score         Freq     Rel Freq    Cum RFreq"
	           "    Perc Rank\n\n");
    for(j=0;j<s->nsv;j++)
      std::fprintf(fp,"%8d%8.3f %12.5f %12.5f %12.5f %12.5f\n",
	          j,score(j,s->minv,s->incv),s->fd_v[j],s->density_v[j],
			  s->crfd_v[j],s->prd_v[j]);
  }

  std::fprintf(fp,"\n\n");
  for(i=1;i<=70;i++)  std::fprintf(fp,"*");

  return;
} 

/*************************************************************************/

void design_matrix(int nsu, double minu, double incu, 
                   int nsv, double minv, double incv,  
                   int cu, int cv, int cuv, int cpm[][2], int scale, 
                   double **B_raw, double **B)
/*
  Create design matrix.  Code and variable names are for a
  bivariate u*v distribution, where u is rows and v is columns.
  Note that x = u + v = total score, with v = common-item score 
  and u = non-common-item score.  We never create a design matrix
  with x and v where v is internal to x, because there would be
  a great deal of collinearity. 
  
  Input
    nsu = number of score categories for u
    minu = minimum score for u
    incu = increment for u
    nsv = number of score categories for v
    minv = minimum score for v
    incv = increment for v
    cu = number of degrees of smoothing for u
    cv = number of degrees of smoothing for v
    cuv = number of cross-product moments
    cpm[cuv-1][2] = zero-offset matrix designating cross-
                    product moments.  Example: let cuv = 3,
                    and the desired cross-product moments be 
                    (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                    Then cpm[0] = 1,1; cpm[1] = 1,2; and
                    cpm[2] = 2,1.  
    scale: 0 --> no scaling; 
           1 --> scale such that each column of B has
                 sum (elements) = 0 and sum (elements^2) = 1
                 
  Output
    B_raw[][] = zero-offset design matrix for raw scores; 
                #rows = (nsu*nsv) and #cols = (cu+cv+cuv);
                space already allocated for B_raw 
    B[][]     = zero-offset design matrix used for solution;
                involves scaling if scale==1 
                #rows = (nsu*nsv) and #cols = (cu+cv+cuv);
                space already allocated for B 
    
  NOTE: For univariate smoothing, set 
        nsv=0, cv=0, cuv=0, cpm = NULL.  In this case, 
        obviously, u plays a generic role 
        (i.e., any single variable such as x or y)
        
  NOTE: Usually for bivariate smoothing in equating, set cuv=1
        and cpm[0] = 1,1. Otherwise, the conditional 
        distributions are not necessarily stochastically ordered
        (see Rosenbaum & Thayer, 1987, p. 46).

  NOTE: Using minu!=0 and incu!=1 changes both the B and 
        B_raw matrices; similarly for minv and incv.  
        If convergence not achieved, it may be wise to set
        minu = 0 and incu = 1 (same for minv and incv)

  Function calls other than C or NR utilities:  
    score()
	runerror()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,j,k,cell,
      ncols = cu+cv+cuv,           /* # columns in design matrix */
      ncells = (nsv>0) ? nsu*nsv : nsu;  /* # rows in design mat */
      
  double mn=0.,                        /* mean for a column of B */
         ss=0.;   /* sum of squared deviations for a column of B */

  double iscore,                           /* score(i,minu,incu) */
         jscore;                           /* score(j,minv,incv) */
  
  /* univariate LL smoothing */
  
  if(nsv==0){    
    for(i=0;i<nsu;i++){
      iscore = score(i,minu,incu);
      for(k=1;k<=cu;k++){ 
        B_raw[i][k-1] = std::pow(iscore, static_cast<double>(k));
        B[i][k-1] = B_raw[i][k-1];
      }
    }
    goto SCALING;
  }
  
  if(cv==0)
    runerror("cv or nsv misspecified");    
  
  /* bivariate LL smoothing */
  
  cell = 0;                                   /* u polynomials */
  for(i=0;i<nsu;i++){
    iscore = score(i,minu,incu);
    for(j=0;j<nsv;j++){
      for(k=1;k<=cu;k++) 
        B_raw[cell][k-1] = std::pow(iscore, static_cast<double>(k));
      cell++;
    }
  }
                                       
  for(j=0;j<nsv;j++){                          /* v polynomials */
    jscore = score(j,minv,incv);
    for(i=0;i<nsu;i++)
      for(k=1;k<=cv;k++) 
        B_raw[j + i*nsv][cu+k-1] = std::pow(jscore, static_cast<double>(k));
  }

  if(cuv>0 && cpm==nullptr)
    runerror("cuv of cpm misspecified");
    
  cell = 0;                       /* cross-product polynomials */
  for(i=0;i<nsu;i++)
    for(j=0;j<nsv;j++){
      for(k=0;k<cuv;k++) 
        B_raw[cell][cu+cv+k] = B_raw[cell][cpm[k][0] - 1]*
                               B_raw[cell][cu + cpm[k][1] - 1];
      cell++;
    }

  
  for(i=0;i<ncells;i++)                     /* copy B_raw to B */
    for(j=0;j<ncols;j++)    
      B[i][j] = B_raw[i][j];

  /* scaling such that for each column of B,
     sum (elements) = 0 and sum (elements^2) = 1 */
  
  SCALING:
  
  if(scale==0) return;                           /* no scaling */
  
  for(j=0;j<ncols;j++){
    mn = ss = 0.;
    for(i=0;i<ncells;i++){
      mn += B[i][j];
      ss += B[i][j]*B[i][j];      
    }
    mn /= ncells;
    ss -= ncells*mn*mn;
      
    for(i=0;i<ncells;i++) B[i][j] = (B[i][j] - mn)/std::sqrt(ss); 
  }
  
  return;

}

/***************************************************************/

void get_nct_bfd(int anchor, int nsx, int nsv, double **bfd, 
                 double *nct)
/*
  Get nct[] from bfd[][]

  Convert bivariate fd (xv->bfd[][]) to a vector nct[] with row j
  elements followed by row j+1 elements.  Note that for an 
  internal anchor the rows are x = u + v and the cols are v,
  which means that there are structural zeros, and we want nct[]
  to contain only the u*v elements.  See example below in which
  - indicates a structural 0 and the within-matrix numbers are 
  the cell locations in nct[].

                 v                                   v
            0  1  2  3                          0  1  2  3
       ---------------                      --------------
       0 |  0  -  -  -                      0|  0  1  2  3
       1 |  4  1  -  -    nsx = 9           1|  4  5  6  7
       2 |  8  5  2  -    nsv = 4 -->     u 2|  8  9 10 11  
       3 | 12  9  6  3    nsu = 6           3| 12 13 14 15
     x 4 | 16 13 10  7                      4| 16 17 18 19
       5 | 20 17 14 11                      5| 20 21 22 23
       6 |  - 21 18 15
       7 |  -  - 22 19
       8 |  -  -  - 23
 

  Input
    anchor : 0 --> external; 1 --> internal
    nsx = number of score categories for total scores (x)
    nsv = number of score categories for common-item scores (v)
    bfd[][] = bivariate freq dist for x and v

  Output
    nct[] = vector version of bfd[][], where nct[] is "collaped",
            as discussed above, if anchor is internal.
            Assumes space already allocated for nct[]
 
  Function calls other than C or NR utilities: None 

  R. L. Brennan

  Date of last revision: 6/30/08      
*/
{
  int i,j,cells,
      nsu; /* number of score categories for non-common items */

  if(!anchor){                             /* external anchor */
    cells = 0;
    for(i=0;i<nsx;i++)                              /* rows x */
      for(j=0;j<nsv;j++)                            /* cols v */
        nct[cells++] = bfd[i][j]; 
  }
  else{                                    /* internal anchor */
    nsu = nsx - nsv + 1;
    for(i=0;i<nsx;i++)                              /* rows x */
      for(j=0;j<nsv;j++){                           /* cols v */
        if(i<j || i>nsu-1+j) continue;
        else nct[(i-j)*nsv+j] = bfd[i][j];    
      }
  }

  return;
}
 
/****************************************************************/

void get_bfd_mct(int anchor, int nsx, int nsv, double *mct, 
                 double **bfd)
/*
  Get bfd[][] from mct[]

  For an external anchor, directly convert the row major 
  vector mct[nsu*nsv] to bfd[nsu][nsv]. For an 
  internal anchor, convert mct[nsu*nsv] to
  bfd[nsu+nsv][nsv] be adding structural zeros.
  See comments for get_nct_bfd().

  Input
    anchor : 0 --> external; 1 --> internal
    nsx = number of score categories for total scores (x)
    nsv = number of score categories for common-item scores (v)
    mct[] = row-major vector for fitted bivariate frequencies
            for non-common items by common items

  Output
    bfd[][] = fitted bivariate frequencies for total scores by
              common-item scores. That is mct[] is mapped into 
              bfd[][] with structural zeros added if anchor is 
              internal (see comments for get_nct_bfd()). 
              Assumes space already allocated for bfd[][]
 
  Function calls other than C or NR utilities: None 

  R. L. Brennan

  Date of last revision: 6/30/08 
*/
{
  int i,j,cells,
      nsu; /* number of score categories for non-common items */

  if(!anchor){                             /* external anchor */
    cells = 0;
    for(i=0;i<nsx;i++)                              /* rows x */
      for(j=0;j<nsv;j++)                            /* cols v */ 
        bfd[i][j] = mct[cells++];
  }
  else{                                    /* internal anchor */
    nsu = nsx - nsv + 1;
    for(i=0;i<nsx;i++)                             /* rows x */
      for(j=0;j<nsv;j++){                           /* cols v */
        if(i<j || i>nsu-1+j) bfd[i][j] = 0.;
        else bfd[i][j] = mct[(i-j)*nsv+j];     
      }
  }
  
  return;
}

/****************************************************************/

void get_BtSmB(double **B, double *m, int ns, int nc, double N,
               double **BtSmB)
/*
  Input:
    B[][] = design matrix (ns x nc)
    m[]   = fitted frequencies (ns x 1)
    ns    = number of score categories (rows in design matrix)
    nc    = number of columns in design matrix
    N     = total of all frequencies

  Output:
    BtSmB[][] = Bt x Sm x B (nc x nc)
              = minus the 2nd derivative of log-likelihood
               (Eq. 22 and 32 in Holland & Thayer, 1987)

  Function calls other than C or NR utilities: None 

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j,t;
  double *Bm;

  Bm = dvector(0,nc-1);

  for(i=0;i<nc;i++){ 
    Bm[i] = 0.;
    for(t=0;t<ns;t++) Bm[i] += B[t][i]*m[t];  
  }

  for(i=0;i<nc;i++)
    for(j=i;j<nc;j++){
      BtSmB[i][j] = 0.;
      for(t=0;t<ns;t++) BtSmB[i][j] += B[t][i]*B[t][j]*m[t];
      BtSmB[i][j] -= Bm[i]*Bm[j]/N;
      BtSmB[j][i] = BtSmB[i][j];
    }

  free_dvector(Bm,0,nc-1);  
}

/****************************************************************/

void get_Btnm(double **B, double *n, double *m, int ns, int nc, 
              double *Btnm)
/*
  Input:
    B[][] = design matrix (ns x nc)
    n[]   = actual frequencies (ns x 1)
    m[]   = fitted frequencies (ns x 1)
    ns    = number of score categories (rows in design matrix)
    nc    = number of columns in design matrix

  Output:
    Btnm[] = 1st derivative of log-likelihood (nc x 1)
             (Eq. 19 and 33 in Holland & Thayer, 1987)

  Function calls other than C or NR utilities: None 

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,t;

  for(i=0;i<nc;i++){ 
    Btnm[i] = 0.;
    for(t=0;t<ns;t++) Btnm[i] += B[t][i]*(n[t] - m[t]);  
  }

}

/****************************************************************/

void get_Beta0(double **B, double *n, double N, int ns, int nc,
               double *Beta0, FILE *fp)
/*
  Input:
    B[][] = zero-offset design matrix (ns x nc) 
    n[]   = zero-offset actual frequencies (ns x 1)
    N     = total of frequencies
    ns    = number of score categories (rows in design matrix)
    nc    = number of columns in design matrix
    fp    = output file pointer for debugging
            (NULL --> no output)

  Output:
    Beta0[] = zero-offset initial values of Beta (nc x 1);
              space already allocated;
              based on Eq 49 (and next line) of Holland and
                Thayer (1987) -- abbreviated H&T below

  Function calls other than C or NR utilities: 
    Print_vector()
    Print_matrix() 

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i,j;
  double *a,                               /* see third line of p. 15 in  H&T */
         **BtSaB,                   /*first term on left side of Eq 49 in H&T */
         *BtSaloga;                                         /* Bt x Sa x loga */

  double Baloga, Ba, aloga=0;

  double **left,     /* left[1...nc][1...nc] set to BtSaB[0...nc-1][0...nc-1] */
	     *right,                    /* right[1..nc] set to BtSaloga[0...nc-1] */
		 d;                                           /* +1 or -1 in ludcmp() */
  int *indx;                            /* row permutation vector in ludcmp() */

  a = dvector(0,ns-1);
  BtSaB = dmatrix(0,nc-1,0,nc-1);
  BtSaloga = dvector(0,nc-1);
 
  /* get a; .8 can be changed to any value in (0,1) */

  for(i=0;i<ns;i++) a[i] = .8*n[i] + .2*N/ns;

       if(fp) 
         Print_vector(fp,"a debug",a,ns," "," ");

 /* get BtSaB -- first term on left side of Eq 49 in H&T */ 

  get_BtSmB(B, a, ns, nc, N, BtSaB); 
      
       if(fp)
         Print_matrix(fp,"BtSaB debug",BtSaB,nc,nc," "," ");

 /* get right side of Eq 49 in H&T, which is computed using
     Equation 38 in Holland and Thayer (2000) with all mu = 0 */
  
  for(i=0;i<ns;i++) aloga += a[i]*std::log(a[i]);
  for(j=0;j<nc;j++){
    Baloga=0.; Ba=0.; 
    for(i=0;i<ns;i++){
      Baloga += B[i][j]*a[i]*std::log(a[i]);
	  Ba += B[i][j]*a[i]; 
    }
	BtSaloga[j] = Baloga - Ba*aloga/N;
  }

  /* get Beta0 using NR ludcmp() and lubksb() */
  
  left = dmatrix(1,nc,1,nc);  
  right = dvector(1,nc);
  indx = ivector(1,nc);

  for(j=0;j<nc;j++){
    right[j+1] = BtSaloga[j];             /* right[] is one-offset B in Ax=B */
    for(i=0;i<nc;i++)
      left[i+1][j+1] = BtSaB[i][j];      /* left[][] is one-offset A in Ax=B */
  }
   
  er_ludcmp(left,nc,indx,&d);    /* left[][] is LU decomp after function call */  
  er_lubksb(left,nc,indx,right);   /* right[] is solution after function call */
  
  for(j=0;j<nc;j++) Beta0[j] = right[j+1];  /* Beta0 is zero-offset solution */ 

  /* deallocate */
  free_dvector(a,0,ns-1);
  free_dmatrix(BtSaB,0,nc-1,0,nc-1);
  free_dvector(BtSaloga,0,nc-1);
  free_dmatrix(left,1,nc,1,nc);
  free_dvector(right,1,nc);
  free_ivector(indx,1,nc);
}

/***************************************************************/

double get_mct(double **B, double *Beta, double *uin, double N, 
               int ns, int nc, double *m, FILE *fp)
/*
  Get m using Equation 8 in Holland and Thayer (1987)

  Input
    B[][]  = design matrix (ns x nc)
    Beta[] = parameter estimates (nc x 1)
    uin[]  = u constants (ns x 1); if NULL, set elements to 0 
    N      = total of frequencies
    ns     = number of rows of design matrix
    nc     = number of columns of design matrix
    fp     = output file pointer for debugging
             (NULL --> no output)

  Output
    m[]    = estimated frequencies (space already allocated)

  NOTE: DBL_MIN is defined in <cfloat>.  It is the minimum normalized 
	    floating point number. For Visual Studio DBL_MIN = 2.225074 E-308
		log(base e) of DBL_MIN in Visual Studio is -708.3964

  Return ap = alpha' --- see top of p. 3 in H&T

  Function calls other than C or NR utilities: 
    Print_vector() 

  R. L. Brennan

  Date of last revision: 6/30/08
*/
{
  int i;
  double *u=nullptr,                                        /* constant */
         *BBeta,                                   /* B x Beta (ns x 1) */
         ap=0.,                                               /* alpha' */
		 ldmin = std::log(DBL_MIN);/* log of smallest number; see NOTE above */

  BBeta = dvector(0,ns-1);

  if(uin==nullptr){
    u = dvector(0,ns-1);
    for(i=0;i<ns;i++) u[i] = 0.;
  }
  else
    u = uin;

  mmult_a_v(B, Beta, ns, nc, BBeta); 
   
       if(fp)
         Print_vector(fp,"BBeta debug", BBeta,ns," "," ");

  for(i=0;i<ns;i++){ 
	ap += (u[i] + BBeta[i] < ldmin) ? 0. : std::exp(u[i] + BBeta[i]);
  }
  ap = std::log(N) - ((ap < ldmin) ? 0. : std::log(ap));

       if(fp)
         std::fprintf(fp,"\n\nap = %12.5f",ap);

  for(i=0;i<ns;i++) 
    m[i] = (ap + u[i] + BBeta[i] < ldmin) ? 0. : std::exp(ap + u[i] + BBeta[i]);
  
/* deallocate */

  if(!uin) free_dvector(u,0,ns-1);
  free_dvector(BBeta,0,ns-1);

  return ap;                        /* ap is the normalizing constant */
}

/**********************************************************************/

int iteration(FILE *fp, double **B, double **B_raw, double *nct, double N, 
              double *uin, int ns, int cu, int cv, int cuv, int cpm[][2], 
              int max_nit, int ctype, int Btype, double crit,
              double *Beta, double *mct, double *n_mts, double *m_mts, 
              double *n_mts_raw, double *m_mts_raw, double *lrc, int *nzero, double *ap)
/*
  Iterate to a solution for log-linear models under a mutinomial 
  distribution as discussed in Holland and Thayer (1987), abbreviated
  H&T here.

  This function is called from Smooth_ULL() and Smooth_BLL()

  Input
    fp      = file pointer for output
              (if fp==NULL) no output written
    B[][]   = design matrix (ns x nc) where nc = cu + cv + cuv
              used for solution
    B_raw[][]   = design matrix for raw scores;
                  used to get central moments               
    nct[]   = frequencies (ns x 1)
    N       = total of frequencies
    uin[]   = constant vector (see Eq 8 in H&T);
              if NULL, all elements set to 0
    ns      = number of frequencies (rows in design matrix)
    cu      = number of degrees of smoothing for u
    cv      = number of degrees of smoothing for v
    cuv     = number of cross-product moments
    cpm[cuv-1][2] = zero-offset matrix designating cross-
                    product moments.  Example: let cuv = 3,
                    and the desired cross-product moments be 
                    (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                    Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                    cpm[2] = {2,1}. 
    max_nit = maximum number of iterations
    ctype = comparison type for criterion:
            0 --> absolute; 1 --> relative
    Btype = type of moments for criterion mathching:
            0 --> use B (could be scaled or unscaled as indicated
                  in design_matrix()) --- see note below
            1 --> use B_raw and central moments based on it
    crit = criterion.  See crit_mts() for discussion 
    
  Output
    Beta[]  = coefficients for variables (columns) of B (nc x 1)
    mct[]   = fitted frequencies (ns x 1) 
    n_mts[] = actual moments based on B (nc x 1)
    m_mts[] = fitted moments based on B (nc x 1)
    n_mts_raw[] = actual central moments based on B_raw (nc x 1)
    m_mts_raw[] = fitted central moments based on B_raw (nc x 1)  
    *lrc        = likelihood-ratio chi-square; see Agresti, 2007, p. 37 
    *nzero      = number of times that nct[i]/mct[i] involves 
                  at least one zero frequency; df correction for *lrc  
	*ap     = alpha (used as normalizing constant in CLL---added by TW)

    NOTE: space already allocated for 
          Beta[], mct[], n_mts[], m_mts[], n_mts_raw, m_mts_raw;
          lrc and nzero declared in calling function

    NOTE: if scale==0 for design_matrix(),
          then n_mts and m_mts are non-central moments;
          if scale==1 for design_matrix(),
          then n_mts and m_mts are analogous to central moments

  Return: number of iterations to convergence 

  Function calls other than C or NR utilities: 
    get_LLmoments()
    Print_iteration_heading() 
    get_Beta0()
    get_mct()
    crit_mts()
    get_BtSmB()
    er_ludcmp()
	er_lubksb()
    get_Btnm()
    mmult_a_v()
    runerror()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{ 
  int i,j,
      nc = cu + cv + cuv,              /* number of columns of design matrix */
      nit,                                               /* iteration number */
      pcells = (ns>250) ? 250: ns;           /* print no more than 250 cells */
                  
  double **BtSmB,                 /* Bt x Sm x B --- see Eq 30 and 32 in H&T */
         *Btnm,                  /* Bt x (n - m) --- see Eq 30 and 33 in H&T */
         *delta,                                         /* see Eq 31 in H&T */
         zero = 0.;

  double **left,     /* left[1...nc][1...nc] set to BtSmB[0...nc-1][0...nc-1] */
	     *right,                        /* right[1..nc] set to Btnm[0...nc-1] */
		 d;                                           /* +1 or -1 in ludcmp() */
  int *indx;                            /* row permutation vector in ludcmp() */

  BtSmB = dmatrix(0,nc-1,0,nc-1);
  Btnm = dvector(0,nc-1);
  delta = dvector(0,nc-1);

  get_LLmoments(B, B_raw, nct, N, ns, cu, cv, cuv, cpm, 
                n_mts, n_mts_raw);                       /* actual mts */

  if(fp) Print_iteration_heading(fp, ns, nc, nct, n_mts, n_mts_raw,
                                 ctype, Btype, crit); 
  get_Beta0(B, nct, N, ns, nc, Beta, nullptr);            /* initial Beta */ 
 

  /***** begin iteration loop *****/ 
  for(nit=0;nit<=max_nit;nit++){ 
    *lrc = 0.;  *nzero = 0;       
	 
    *ap = get_mct(B, Beta, uin, N, ns, nc, mct, nullptr); /* get ap & mct */
	/* ap is the normalizing constant which is needed in CLL */

    for(i=0;i<ns;i++)         
      if(nct[i]!=0. && mct[i]!=0.) *lrc += nct[i]*std::log(nct[i]/mct[i]);
      else (*nzero)++;         
    *lrc *= 2;                          /* likelihood-ratio chi-square */

    get_LLmoments(B, B_raw, mct, N, ns, cu, cv, cuv, cpm,
                  m_mts, m_mts_raw);                     /* fitted mts */

    if(fp){                  /* print fitted-results for iteration nit */
      std::fprintf(fp,"\n        %2d    ",nit);
      std::fprintf(fp,"  %10.5f",*ap);
      std::fprintf(fp,"  %10.5f",(uin==nullptr) ? zero : uin[i]);
      for(i=0;i<nc;i++) std::fprintf(fp,"%12.3f",Beta[i]);
      for(i=0;i<nc;i++) std::fprintf(fp,"%12.5f",m_mts[i]);
      for(i=0;i<nc;i++) std::fprintf(fp,"%12.5f",m_mts_raw[i]);
      std::fprintf(fp,"%12.5f",*lrc);
      std::fprintf(fp,"%7d",*nzero);
      for(i=0;i<pcells;i++) std::fprintf(fp,"%12.5f",mct[i]);
      if(ns>250) std::fprintf(fp,"  ...");
    }
      
    if(crit_mts(nc, cu, ctype, Btype, 
                  (Btype==0) ? n_mts : n_mts_raw, 
                  (Btype==0) ? m_mts : m_mts_raw,
                  crit)) break;                             /* mts criterion */

    get_BtSmB(B, mct, ns, nc, N, BtSmB);                        /* get BtSmB */ 
    get_Btnm(B, nct, mct, ns, nc, Btnm);                         /* get Btnm */

	/* solve for delta^r in ER Equation 10.14 == H&T Equation 30*/

    left = dmatrix(1,nc,1,nc);  
    right = dvector(1,nc);
    indx = ivector(1,nc);

    for(j=0;j<nc;j++){
      right[j+1] = Btnm[j];                /* right[] is one-offset B in Ax=B */
      for(i=0;i<nc;i++)
        left[i+1][j+1] = BtSmB[i][j];      /* left[][] is one-offset A in Ax=B */
    } 
    
    er_ludcmp(left,nc,indx,&d);   /* left[][] is LU decomp after function call */
    er_lubksb(left,nc,indx,right);  /* right[] is solution after function call */
    
	/* get delta[], which is zero-offset version of right[] */ 
	for(j=0;j<nc;j++) delta[j] = right[j+1]; 

	/* get new Beta */ 
    for(j=0;j<nc;j++) Beta[j] += delta[j];                

  } 

  /***** end iteration loop *****/                                             

  if(nit>max_nit)
    runerror("Criterion not satisfied after max_nit iterations");                                            

  free_dmatrix(BtSmB,0,nc-1,0,nc-1);
  free_dvector(Btnm,0,nc-1);
  free_dvector(delta,0,nc-1);
  free_dmatrix(left,1,nc,1,nc);
  free_dvector(right,1,nc);
  free_ivector(indx,1,nc);

  return nit;
}

/**********************************************************************/

void get_LLmoments(double **B, double **B_raw, double *f, double N, 
                   int ns, int cu, int cv, int cuv, int cpm[][2], 
                   double *mts, double *mts_raw)
/*
  Get moments based on B and B_raw design matrices for log-linear model

  Input
    B[][] = design matrix (ns x nc)
            could be based on scale = 0 or scale = 1 as specified
            in design_matrix()
    B_raw[][] = design matrix (ns x nc) for raw scores
                (e.g., first column is min, min+inc, ...)
    f[]   = frequencies (ns x 1) (could be nct[] or mct[])
    ns    = number of frequencies (rows in design matrix)
    cu    = number of degrees of smoothing for u
    cv    = number of degrees of smoothing for v
    cuv   = number of cross-product moments
    cpm[cuv-1][2] = zero-offset matrix designating cross-
                    product moments.  Example: let cuv = 3,
                    and the desired cross-product moments be 
                    (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                    Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                    cpm[2] = {2,1}.  
  Output
    mts[] = moments (nc x 1) based on B (space already allocated)
    mts_raw[] = central moments (nc x 1) based on B_raw 
                (space already allocated);
                e.g., for j = 2,...(cu-1)), the (j+1)th central moment is 
                      mts_raw[j] = 
                      sum_i(B_raw[i][j] - xbar)^(j+1)*f[i]/sd^((j+1)),
                      where xbar = sum_i(B_raw[i][0]*f[i]).
                      Note that mts_raw[0] = 0 by construction, 
                      and we replace it with xbar 
                      (i.e., we set mts_raw[0] = xbar);
                      simlarly, we set mts_raw[1] = sd

  Function calls other than C or NR utilities: none

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i,j,
      nc = cu + cv + cuv;       /* number of columns in design matrix */
  double mnu, sdu, mnv, sdv, 
         *rf;                                   /* relative frequency */

  rf = dvector(0,ns-1);
  for(i=0;i<ns;i++) rf[i] = f[i]/N;

  for(j=0;j<nc;j++) mts[j] = mts_raw[j] = 0.;       /* initialization */

  /*** moments based on B ***/
                                        
  for(j=0;j<nc;j++)
    for(i=0;i<ns;i++)
      mts[j] += B[i][j]*rf[i];

  /*** "typical" central moments based on B_raw ***/

       /* for cu columns associated with u;  scores are in column 0 */ 

  mnu = sdu = 0.;
  for(i=0;i<ns;i++){
    mnu += B_raw[i][0]*rf[i];                  
    sdu += B_raw[i][0]*B_raw[i][0]*rf[i];     
  }
  sdu = std::sqrt(sdu - mnu*mnu);

  if(cu>=1) mts_raw[0] = mnu;
  if(cu>=2) mts_raw[1] = sdu;
  for(j=2;j<cu;j++){                             /* start with skew */ 
    for(i=0;i<ns;i++) 
      mts_raw[j] += std::pow(B_raw[i][0] - mnu,j+1)*rf[i];
    mts_raw[j] /= std::pow(sdu,static_cast<double>(j+1));
  }

        /* for cv columns associated with v; scores are in column cu*/ 

  if(cv>0){

    mnv = sdv = 0.;
    for(i=0;i<ns;i++){
      mnv += B_raw[i][cu]*rf[i];                 
      sdv += B_raw[i][cu]*B_raw[i][cu]*rf[i];     
    }
    sdv = std::sqrt(sdv - mnv*mnv);

    if(cv>=1) mts_raw[cu] = mnv;
    if(cv>=2) mts_raw[cu+1] = sdv;
    for(j=cu+2;j<cu+cv;j++){                     /* start with skew */
      for(i=0;i<ns;i++) 
        mts_raw[j] += std::pow(B_raw[i][cu] - mnv,j-cu+1)*rf[i];
      mts_raw[j] /= std::pow(sdv,static_cast<double>(j-cu+1));
    }

  }

       /* for cuv columns associated with cross products;
          scores are in columns 0 and cu */ 

  if(cuv>0 && cu>=2 && cv>=2){

    for(j=cu+cv;j<nc;j++){
      for(i=0;i<ns;i++)
        mts_raw[j] += std::pow(B_raw[i][0]  - mnu,cpm[j-cu-cv][0])*
                      std::pow(B_raw[i][cu] - mnv,cpm[j-cu-cv][1])*rf[i];
      mts_raw[j] /= std::pow(sdu,cpm[j-cu-cv][0])*std::pow(sdv,cpm[j-cu-cv][1]);
    }

  }

  free_dvector(rf,0,ns-1);

  return;  

}

/**********************************************************************/

int crit_mts(int nc, int cu, int ctype, int Btype,
             double *n_mts, double *m_mts, double crit)
/*
  Moments criterion for iteration loop

  Input
    nc = number of cmoments
    cu = number of moments for u
    ctype = comparison type for criterion:
            0 --> absolute; 1 --> relative
    Btype = type of moments for criterion mathching:
            0 --> moments based on B 
			      (if scale = 0, design matrix is based on raw scores,
			       which means that the moments are based on raw scores;
			       if scale = 1, design matrix is based on scaled raw scores,
			       which means that the moments are based on scaled raw scores) 
            1 -->  moments based on B_raw, whether scale is 0 or 1
			NOTE: see design_matrix() for comments about scale
    n_mts[] = moments based on actual frequencies
    m_mts[] = moments based on fitted frequencies
    crit = criterion  (see notes below)

  return: 0 --> criterion not met ---- continue iteration
          1 --> criterion met --- end iteration

  NOTES.

  If ctype==0 (absolute comparison), criterion is met if
  |n_mts[i] - m_mts[i]| <= crit for all nc moments.

  If ctype = 1 (relative comparison), criterion is met if
  |(n_mts[i] - m_mts[i])/n_mts[i]| <= crit for all nc moments.  
  An error occurs if n_mts[i] = 0.

  If Btype==1, the first central moment for both u and v
  is 0 whether the frequencies are actual (n) or fitted (m).
  These are the moments associated with column 0 and cu in B[][].
  Further, when Btype==1, these moments are  replaced by u and v 
  means, respectively, in get LLmoments().Therefore, when Btype==1,
  no test is made for |n_mts[0] - m_mts[0]| or for
  |n_mts[cu] - m_mts[cu]|. 

  Function calls other than C or NR utilities:
    runerror()

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
  int i;
  
  if(Btype==0){            /* Btype == 0 --> based on B matrix */
    for(i=0;i<nc;i++){
      if(ctype==0){                     /* absolute comparison */
        if(std::fabs(n_mts[i] - m_mts[i]) > crit) return 0;
      }
      else{                             /* relative comparison */
        if(n_mts[i]==0.) 
          runerror("\n\nRelative criterion and n_mts[] = 0");
        if(std::fabs((n_mts[i] - m_mts[i])/n_mts[i]) > crit) return 0;
      }
    }
  }
  else{ /* Btype==1; based on B_raw matrix and central moments */
    for(i=0;i<nc;i++){
      if(i==0 || i==cu) continue;
      if(ctype==0){                     /* absolute comparison */
        if(std::fabs(n_mts[i] - m_mts[i]) > crit) return 0;
      }
      else{                             /* relative comparison */
        if(n_mts[i]==0.) 
          runerror("\n\nRelative criterion and n_mts[] = 0");
        if(std::fabs((n_mts[i] - m_mts[i])/n_mts[i]) > crit) return 0;
      }
    }
  }

  return 1;                                   /* criterion met */

}

/***************************************************************/ 

void Print_iteration_heading(FILE *fp, int ns, int nc, double *nct,
                             double *n_mts, double *n_mts_raw,
                             int ctype, int Btype, double crit)
/*
  For log-linear iterations print (a) actual results based on nct[]
  and (b) heading for fitted results

  Input
    fp          = file pointer for output
    ns          = number of score categories 
                  (number of rows in design matrix)
    nc          = number of columns in design matrix
    nct[]       = actual frequencies
    n_mts[]     = moments based on nct[] and design matrix B
    n_mts_raw[] = moments based on nct[] and design matrix B_raw
    ctype = comparison type for criterion:
            0 --> absolute; 1 --> relative
    Btype = type of design matrix and, hence, type of moments 
            for criterion:
            0 --> use B (could be scaled or unscaled as indicated
                  in design_matrix()) 
            1 --> use B_raw and central moments based on it
    crit = criterion.  See crit_mts() for discussion 

  Function calls other than C or NR utilities: none

  R. L. Brennan

  Date of last revision: 6/30/08


*/
{
  int i,
      pcells = (ns>250) ? 250: ns; 

    /* print headings */

    std::fprintf(fp,"\n\n\nRESULTS FOR ITERATIONS");

    std::fprintf(fp,"\n\n              ");
    for(i=1;i<=(12*(nc+2));i++) std::fprintf(fp," ");

    for(i=1;i<=((12*nc)-18)/2;i++) std::fprintf(fp," ");
    std::fprintf(fp,"Moments Based on B");
    for(i=1;i<=((12*nc)-18)/2;i++) std::fprintf(fp," ");

    for(i=1;i<=((12*nc)-30)/2;i++) std::fprintf(fp," ");
    std::fprintf(fp,"Central Moments Based on B_raw");
    for(i=1;i<=((12*nc)-30)/2;i++) std::fprintf(fp," ");

    std::fprintf(fp,"\n              ");
    for(i=1;i<=(12*(nc+2));i++) std::fprintf(fp," ");
    std::fprintf(fp," ");
    for(i=1;i<=12*nc-1;i++) std::fprintf(fp,"-");
    std::fprintf(fp," ");
    for(i=1;i<=12*nc-1;i++) std::fprintf(fp,"-");

    /* print actual results based on nct[] */

    std::fprintf(fp,"\n              ");
    for(i=1;i<=(12*(nc+2));i++) std::fprintf(fp," ");
    for(i=0;i<nc;i++) std::fprintf(fp,"       n[%2d]",i+1);
    for(i=0;i<nc;i++) std::fprintf(fp,"   n_raw[%2d]",i+1); 
    for(i=0;i<19;i++) std::fprintf(fp," ");    
    for(i=0;i<pcells;i++) std::fprintf(fp,"    nct[%3d]",i);
    if(ns>250) std::fprintf(fp,"  ...");
    std::fprintf(fp,"\n\n        Actual");
    for(i=1;i<=(12*(nc+2));i++) std::fprintf(fp," ");
    for(i=0;i<nc;i++) std::fprintf(fp,"%12.5f",n_mts[i]);
    for(i=0;i<nc;i++) std::fprintf(fp,"%12.5f",n_mts_raw[i]);
    for(i=0;i<19;i++) std::fprintf(fp," ");
    for(i=0;i<pcells;i++) std::fprintf(fp,"%12.5f",nct[i]);
    if(ns>250) std::fprintf(fp,"  ...");
 
    /* print criterion */

    std::fprintf(fp,"\n\n    Criterion:  ");
    if(Btype==0 && ctype==0)
      std::fprintf(fp,"|(n[i] - m[i]| <= %15.10f",crit); 
    else if(Btype==0 && ctype==1)
      std::fprintf(fp,"|n[i] - m[i])/n[i]| <= %15.10f",crit);
    else if(Btype==1 && ctype==0)
      std::fprintf(fp,"|(n_raw[i] - m_raw[i]| <= %15.10f",crit); 
    else 
      std::fprintf(fp,"|n_raw[i] - m_raw[i])/n_raw[i]| <= %15.10f", crit);
    std::fprintf(fp,"  for i = 1 to %d",nc);

    /* print fitted-results heading for iterations */

    std::fprintf(fp,"\n\n     Iteration");
    std::fprintf(fp,"          ap");
    std::fprintf(fp,"           u");
    for(i=0;i<nc;i++) std::fprintf(fp,"    Beta[%2d]",i+1);
    for(i=0;i<nc;i++) std::fprintf(fp,"       m[%2d]",i+1);
    for(i=0;i<nc;i++) std::fprintf(fp,"   m_raw[%2d]",i+1);
    std::fprintf(fp,"     LR Chi2");
    std::fprintf(fp,"  nzero");
    for(i=0;i<pcells;i++) std::fprintf(fp,"    mct[%3d]",i);
    if(ns>250) std::fprintf(fp,"  ...");
    std::fprintf(fp,"\n");
}

/*********************************************************************/

void mtranspose(double **a, int nr, int nc, double **at)
/* Transpose a matrix.
 
   Input
     a = matrix (nr x nc)
     nr = number of rows in a
     nc = number of cols in a    
     
   Output
     at = transposed matrix (nc x nr)

  Function calls other than C or NR utilities: none

  R. L. Brennan

  Date of last revision: 6/30/08

*/
{
   int i,j;
   
   for(i=0; i<nr; i++)                             /* row of a */
     for(j=0; j<nc; j++)                           /* col of a */
       at[j][i] = a[i][j];
                                
}

/***************************************************************/

void mmult_a_b(double **a, double **b, 
               int nra, int nca, int nrb, int ncb,
               double **s)
/* Multiply two matrices (a and b), and store results 
   in s = a * b;  Error if nca != nrb
   
   Input
     a = matrix (nra x nca)
     b = matrix (nrb x ncb) 
     nra = number of rows for a
     nca = number of cols for a
     nrb = number of rows for b
     ncb = number of cols for b 
     
   Output
     s = matrix (space already allocated --- see note, below)

     NOTE:  s = a * b is (nra * ncb) with
            rows from 0...(nra-1) and
            cols from 0...(ncb-1).  

  Function calls other than C or NR utilities: none

  R. L. Brennan

  Date of last revision: 6/30/08
                        
*/
{
  int i,j,k;
   
  if(nca!=nrb)
     runerror("\n\nError in mmult_a_b()");
                                        
  for(i=0; i<nra; i++)                   /* row of a; row of s */
    for(j=0; j<ncb; j++){                /* col of b; col of s */
      s[i][j] = 0.;                      
      for(k=0; k<nca; k++)                      
          s[i][j] += a[i][k]*b[k][j];  /* across row; down col */            
      }
           
  return;
                                
}
         
/***************************************************************/

void mmult_a_v(double **a, double *v, int nr, int nc, double *s) 
/* Multiply matrix a and vector v;
   assumes a is nr x nc, b has nc elements, s has nr elements;
   if we view v as a column vector, then 
   s = a x v is a column vector with nr elements
   
   Input
     a = matrix (nr x nc)
     v = (column) vector (nc elements)
     nr = number of rows in a 
     nc = number of cols in a   
     
   Output
     s = (column) vector (nr elements) 

  Function calls other than C or NR utilities: none

  R. L. Brennan

  Date of last revision: 6/30/08
  
*/
{
  int i,k;
   
  for(k=0; k<nr; k++){                     /* elements of s */                   
    s[k] = 0.;                            
    for(i=0; i<nc; i++)               
      s[k] += a[k][i]*v[i];
    }
    
  return;
}

/***************************************************************/