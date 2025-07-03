/*
  CG_NoSmooth.cpp 

  Wrapper, print, and linear functions for common-item non-equivalent groups 
  design with no smoothing
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
#include "CG_EquiEquate.h"
#include "CG_NoSmooth.h"
#include <cstring>
#include <cstdio>
#include <cmath>

/*******************************************************************************/

void Wrapper_CN(char design, char method, char smoothing, 
               double w1, int anchor, double rv1, double rv2,
               struct BSTATS *xv, struct BSTATS *yv, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for getting mean, linear, or equipercentile equating for CINEG design
  with no smoothing. Equipercentile equating includes frequency estimation with 
  Braun-Holland (linear) results, modified frequency estimation with 
  Braun-Holland (linear) results, and chained equipercentile equating

  Purposes of Wrapper_CN() are to:
      (i) allocate space for r->eraw[][] and r->mts[0][] 
     (ii) populate elements of inall, including the 
          determination of proportional weights, if requested
    (iii) get equating results and store them in r
     (iv) get moments
    
  Assumes that in xv, score 1 is for x and score 2 is for v
  Assumes that in yv, score 1 is for y and score 2 is for v
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep
  
  Input:
  
    design = 'C' (CINEG)
    method:  'M' = mean, 
             'L' = linear (Tucker, Levine-obs, Levine-true, chained)
             'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
             'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
             'G' = FE + BH-FE + MFE + BH-MFE
             'C' = Chained
			 'H' = FE + BH-FE + Chained
             'A' = FE + BH-FE + MFE + BH-MFE + Chained
              
    smoothing = 'N' (none)  
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
    rep = replication number for bootstrap; should be set to 0
          for actual equating; 
   
    NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
    
  Output:
    
    struct PDATA inall:   populates selected values of inall 
    
    struct ERAW_RESULTS r          

      msx[] = means for x for synthetic pop 
      msy[] = means for y for synthetic pop 
      dssx[] = sd's for x for synthetic pop 
      ssy[] = sd's for y for synthetic pop 
      gamma1[] = gamma's for pop 1 
      gamma2[] = gamma's for pop 2 
      a[] = slopes
      b[] = intercepts
      eraw[][]:  equated raw scores
      mts[][]:  moments for equated raw scores  
      
      NOTE: vectors in ERAW_RESULTS for nm = 4 linear methods 
            give results for following methods: 0 - Tucker,
            1 - Levine Observed, 2 - Levine True, 3 - Chained  
          
      NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
            of equated scores; memory allocated here;
            eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
            because we are getting equated raw scores for x.
            eraw[][] stored in this "row"  manner so that 
            Equated_ss() can be used directly on 
            eraw[k] where k is the method number  
          
  NOTE: Whenever method differs, there must be different structures
        passed as struct PDATA and struct ERAW_RESULTS 
    
  NOTE: If Wrapper_CN() is called in a bootstrap loop, then in
        the calling function struct ERAW_RESULTS must be different
        from struct ERAW_RESULTS for the actual equating. 
                                            
  Function calls other than C or NR utilities:
    CI_LinEq() 
    FEorMFE_EE()
    Chained_EE()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08   
*/
{ 
  int i;
                       /* method names --- 10 characters; right justified */
  char *names[] ={"    Tucker", "   Lev Obs", "  Lev True", "  ChainedL",
                  "        FE", "     BH-FE", "       MFE", "    BH-MFE",
                  "  ChainedE"};                   
  
  inall->rep = rep;               /* should be set to 0 for actual equating. */
                    /* Counting of replications done in Wrapper_Bootstrap(), 
             which is why this statement cannot be in the if statement below */ 
                    
  /* allocation and assignments for inall
     Note that for every assignment of the form inall->(var) = xv->(var)
	 or inall->(var) = yv->(var) values vary depending on whether xv is for actual equating or
     a bootstrap sample; all other values are the same for the 
     actual equating and a bootstrap sample */
  
  if(inall->rep == 0){   /* no assignment or stor alloc for bootstrap reps */
    std::strcpy(inall->xfname,xv->fname);
    std::strcpy(inall->yfname,yv->fname);
    inall->xv = xv;
    inall->yv = yv;
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
    if(method == 'M' || method == 'L'){            /* method == 'M' or 'L' */
      inall->nm = 4;
      for(i=0;i<=3;i++) std::strcpy(inall->names[i],names[i]);
    }
    else if(method == 'E'){                               /* method == 'E' */
      inall->nm = 2;
      for(i=4;i<=5;i++) std::strcpy(inall->names[i-4],names[i]);
    }
    else if(method == 'F'){                               /* method == 'F' */
      inall->nm = 2;
      for(i=6;i<=7;i++) std::strcpy(inall->names[i-6],names[i]);
    }
    else if(method == 'G'){                               /* method == 'G' */
      inall->nm = 4;
      for(i=4;i<=7;i++) std::strcpy(inall->names[i-4],names[i]);
    }
    else if(method == 'C'){                               /* method == 'C' */
      inall->nm = 1;
      std::strcpy(inall->names[0],names[8]);
    }
	else if(method == 'H'){                               /* method == 'H' */
      inall->nm = 3;
	  for(i=4;i<=5;i++) std::strcpy(inall->names[i-4],names[i]);
      std::strcpy(inall->names[2],names[8]);
    }
    else{                                                 /* method == 'A' */
      inall->nm = 5;
      for(i=4;i<=8;i++) std::strcpy(inall->names[i-4],names[i]);
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
   
/* Mean and Linear equating results: 

   xv and yv variables relative to variables in CI_LinEq() 

   Function  Function Call
    mnx1      xv->mts1[0]
    sdx1      xv->mts1[1]
    mnv1      xv->mts2[0]
    sdv1      xv->mts2[1]
    covxv1    xv->cov

    mny2      yv->mts1[0]
    sdy2      yv->mts1[1]  
    mnv2      yv->mts2[0]
    sdv2      yv->mts2[1]
    covyv2    yv->cov  
    
    In function, 1 and 2 are populations; in xv and yv, 1 and 2 are variables  
*/

  if(method == 'M' || method == 'L')  
    CI_LinEq(xv->mts1[0], xv->mts1[1], xv->mts2[0], xv->mts2[1], xv->cov,
             yv->mts1[0], yv->mts1[1], yv->mts2[0], yv->mts2[1], yv->cov,
             inall->w1, anchor, method, xv->min1, xv->max1, xv->inc1, inall->nm,           
             r->msx, r->msy, r->ssx, r->ssy,
             r->gamma1, r->gamma2, r->a, r->b, r->eraw); 

/* Equipercentile results, including Braun-Holland (BH) linear. 
   Note: For FE syn densities are in fxs[0] and gys[0]
         For MFE syn densities are in fxs[1] and gys[1] 
         For BH under FE, slope in a[0] and intercept in b[0]
         For BH under MFE, slope in a[1] and intercept in b[1]
*/
                                           /* FE + BH-FE in positions 0 and 1*/
  if(method == 'E' || method == 'G' || method == 'A' || method == 'H')     
    FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1, 
               yv->ns1, yv->min1, yv->max1, yv->inc1,
               xv->bp12, yv->bp12, 0, 0, 
               r->fxs[0], r->gys[0], r->eraw[0],
               &r->a[0], &r->b[0], r->eraw[1]); 

  if(method == 'F')                     /* MFE + BH-MFE in positions 0 and 1 */
    FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1, 
               yv->ns1, yv->min1, yv->max1, yv->inc1,
               xv->bp12, yv->bp12, inall->rv1, inall->rv2, 
               r->fxs[1], r->gys[1], r->eraw[0],
               &r->a[1], &r->b[1], r->eraw[1]);    

  if(method == 'G' || method == 'A')    /* MFE + BH-MFE in positions 2 and 3 */
    FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1, 
               yv->ns1, yv->min1, yv->max1, yv->inc1,
               xv->bp12, yv->bp12, inall->rv1, inall->rv2, 
               r->fxs[1], r->gys[1], r->eraw[2],
               &r->a[1], &r->b[1], r->eraw[3]); 

  /* Chained equipercentile method */
  double* ptr = nullptr;
  if(method == 'C')                                 /* Chained in position 0 */
    ptr = r->eraw[0];
  else if(method == 'A') ptr = r->eraw[4];
  else if(method == 'H') ptr = r->eraw[2];

  if(ptr != nullptr)
    Chained_EE(xv->ns1, xv->prd1, xv->min2, xv->max2, xv->inc2, 
               xv->ns2, xv->crfd2, yv->min1, yv->inc1, yv->ns1,   
               yv->crfd1, yv->crfd2, ptr);  
                        
/* get moments */

  for(i=0;i<=inall->nm-1;i++) 
    MomentsFromFD(xv->min1,xv->max1,xv->inc1,r->eraw[i],xv->fd1,r->mts[i]);
  
}

/******************************************************************************/
          
void CI_LinEq(double mnx1, double sdx1, double mnv1, double sdv1, double covxv1,
              double mny2, double sdy2, double mnv2, double sdv2, double covyv2,
              double w1, int anchor, char method, double min, double max,
              double inc, int nm, double *msx ,double *msy, double *ssx, 
              double *ssy, double *gamma1, double *gamma2, double *a, double *b, 
              double **eraw)
/*
 computes results for common-item linear equating
 of x to scale of y for methods:
   0. Tucker
   1. Levine Observed
   2. Levine True
   3. Chained
            
 Kolen and Brennan (2004) notation used here
 
 ... (Input/Output comments omitted for brevity) ...

  Function calls other than C or NR utilities:
    CI_LinObsEq()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08                                 
*/
{
  int i,j;
  
  if(method=='M'){                                           /* mean equating */
    a[0] = a[1] = a[2] = a[3] = 1.;                       /* slopes are unity */          
    for(i=0;i<=3;i++)
      ssx[i] = ssy[i] = 0.;            /* no synthetic sd's for mean equating */
  }

  msx[2] = msy[2] = ssx[2] = ssy[2] = 0.;       /* For Levine true, synthetic
                                                     population is irrelevant */
/* Tucker */

  gamma1[0] = covxv1/(sdv1*sdv1);
  gamma2[0] = covyv2/(sdv2*sdv2);
  CI_LinObsEq(mnx1,sdx1,mnv1,sdv1,mny2,sdy2,mnv2,sdv2,
              w1,method,gamma1[0],gamma2[0],
              &msx[0],&msy[0],&ssx[0],&ssy[0],&a[0],&b[0]);

/* Levine */

  if(anchor != 0) {                                        /* internal anchor */
    gamma1[1] = sdx1*sdx1/covxv1;
    gamma2[1] = sdy2*sdy2/covyv2;
  }
  else {                                                   /* external anchor */
    gamma1[1] = (sdx1*sdx1 + covxv1)/(sdv1*sdv1 + covxv1);
    gamma2[1] = (sdy2*sdy2 + covyv2)/(sdv2*sdv2 + covyv2);
  }
  
  CI_LinObsEq(mnx1,sdx1,mnv1,sdv1,mny2,sdy2,mnv2,sdv2,  
              w1,method,gamma1[1],gamma2[1],
              &msx[1],&msy[1],&ssx[1],&ssy[1],&a[1],&b[1]); /* Levine lin obs */  
  
  gamma1[2] = gamma1[1];
  gamma2[2] = gamma2[1];
  if(method != 'M') a[2] = gamma2[2]/gamma1[2];      /* Levine lin true slope */
  b[2] = (mny2 - a[2]*mnx1) + gamma2[2]*(mnv1-mnv2); /* Levine lin true inter */
 
/* Chained (Note:  Chained true = Levine true) */

  gamma1[3] = sdx1/sdv1;
  gamma2[3] = sdy2/sdv2; 
  
  CI_LinObsEq(mnx1,sdx1,mnv1,sdv1,mny2,sdy2,mnv2,sdv2,
              w1,method,gamma1[3],gamma2[3],
              &msx[3],&msy[3],&ssx[3],&ssy[3],&a[3],&b[3]); 
              
/* NOTE: The synthetic group means and sd's don't really apply to
   chained linear equating in the sense that the slope and intercept are
   invariant with respect to choice of weights.  However, 
   CI_LinObsEq() is a convenient way to get the slope and intercept */
                            
/* get equated raw scores */               
  
  for(i=0;i<=nm-1;i++)
    for(j=0;j<=loc(max,min,inc);j++)
      eraw[i][j] = b[i] + a[i]*score(j,min,inc); 

}

/*******************************************************************************/

void CI_LinObsEq(double mnx1, double sdx1, double mnv1, double sdv1, 
                 double mny2, double sdy2, double mnv2, double sdv2,
                 double w1, char method, double gamma1, double gamma2,
                 double *msx, double *msy, double *ssx, double *ssy,
                 double *a, double *b)
/* For a linear observed score equating method and a CINEG design,
  get synthetic pop means and sd's, and get slope (*a) and intercept (*b)  
... (comments omitted for brevity) ...                                              
*/
{
  double musx,musy,varsx,varsy,w2=1-w1;
  
  musx = mnx1 - w2*gamma1*(mnv1-mnv2);
  musy = mny2 + w1*gamma2*(mnv1-mnv2);
  
  if(method != 'M'){
    varsx = sdx1*sdx1 - w2*gamma1*gamma1*(sdv1*sdv1-sdv2*sdv2)
                      + w1*w2*gamma1*gamma1*(mnv1-mnv2)*(mnv1-mnv2);
    varsy = sdy2*sdy2 + w1*gamma2*gamma2*(sdv1*sdv1-sdv2*sdv2)
                      + w1*w2*gamma2*gamma2*(mnv1-mnv2)*(mnv1-mnv2);
    *a = std::sqrt(varsy/varsx);
    *ssx = std::sqrt(varsx);
    *ssy = std::sqrt(varsy);
  }
  
  *b = musy - (*a)*musx;
  *msx = musx;
  *msy = musy;
}

/*************************************************************************/

void Print_CN(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
/*
  print CN results (CINEG; mean, linear, equipercentile; no smoothing)
 ... (comments omitted for brevity) ...    
*/
  int i,j,
      jlast = (inall->method=='A') ? 75 : 63;
  
  std::fprintf(fp,"\n\n");
  for(j=1;j<=jlast;j++) std::fprintf(fp,"*");
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  if(inall->method == 'M' || inall->method == 'L'){
    std::fprintf(fp,"%s Equating with CINEG Design (%s Anchor; w1 = %7.5f)\n\n",
               (inall->method == 'M') ? "Mean" : "Linear",
               (inall->anchor != 0) ? "Internal" : "External", inall->w1); 
    std::fprintf(fp,"Input file for xv and pop 1: %s\n",inall->xfname);
    std::fprintf(fp,"Input file for yv and pop 2: %s\n\n",inall->yfname);
    std::fprintf(fp,"                  Method 0:   Method 1:   Method 2:   Method 3:\n");
    std::fprintf(fp,"                 %s  %s  %s  %s\n\n",
                inall->names[0],inall->names[1],inall->names[2],inall->names[3]);
    std::fprintf(fp,"          msx  %12.5f%12.5f\n",r->msx[0],r->msx[1]);
    std::fprintf(fp,"          msy  %12.5f%12.5f\n",r->msy[0],r->msy[1]);
          
    if(inall->method != 'M'){        
      std::fprintf(fp,"          ssx  %12.5f%12.5f\n",r->ssx[0],r->ssx[1]);
      std::fprintf(fp,"          ssy  %12.5f%12.5f\n",r->ssy[0],r->ssy[1]);
    }
          
    std::fprintf(fp,"\n       gamma1  %12.5f%12.5f%12.5f%12.5f\n",
            r->gamma1[0],r->gamma1[1],r->gamma1[2],r->gamma1[3]);
    std::fprintf(fp,"       gamma2  %12.5f%12.5f%12.5f%12.5f\n\n",
            r->gamma2[0],r->gamma2[1],r->gamma2[2],r->gamma2[3]);
    std::fprintf(fp,"      inter b  %12.5f%12.5f%12.5f%12.5f\n",
            r->b[0],r->b[1],r->b[2],r->b[3]);
    std::fprintf(fp,"      slope a  %12.5f%12.5f%12.5f%12.5f\n",
            r->a[0],r->a[1],r->a[2],r->a[3]);
  }
  else{                /* print statements for FE, BH-FE, MFE, BH-MFE, Chained */
    std::fprintf(fp,"Equipercentile Equating with CINEG Design"
               " (%s Anchor; w1 = %7.5f)\n\n",
               (inall->anchor != 0) ? "Internal" : "External", inall->w1); 
    std::fprintf(fp,"Input file for xv and pop 1: %s\n",inall->xfname);
    std::fprintf(fp,"Input file for yv and pop 2: %s\n\n",inall->yfname);

    if(inall->method == 'E' || inall->method == 'G' || inall->method == 'A')
      std::fprintf(fp,"Braun-Holland (BH) Under Frequency Estimation (FE):\n"
                 "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
                 r->b[0], r->a[0]);
    if(inall->method == 'F' || inall->method == 'G' || inall->method == 'A')
      std::fprintf(fp,"Braun-Holland (BH) Under Modified Frequency Estimation (MFE):\n"
                 "    Intercept: b = %12.5f;     Slope: a = %12.5f\n\n",
                 r->b[1], r->a[1]);
  }

  for(i=1;i<=jlast;i++) std::fprintf(fp,"-");

  /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (x)  ");
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
/*********************************************************************************/