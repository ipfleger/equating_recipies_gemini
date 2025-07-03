/* CG_EquiEquate.cpp - contains functions for computing equating results for 
          (a) frequency estimation (FE) and Braun-Holland 
              under FE (if requested);
          (b) modified frequency estimation (MFE) (Wang & Brennan, 2006) 
              and Braun-Holland under MFE (if requested); and
          (c) chained equipercentile equating

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
#include <cmath>

/*****************************************************************************/
     
void FEorMFE_EE(double w1, int internal, int nsv, 
                int nsx, double minx, double maxx, 
                int nsy, double miny, double maxy, double inc,
                double **bxvin, double **byvin, double rv1, double rv2,
                double *fxs, double *gys, double *eraw,
                double *a, double *b, double *erawBH)
{
/*
   Computes results for common-item (CI) equipercentile equating for EITHER
   (a) frequency estimation (FE) and Braun-Holland under FE (if requested); OR
   (b) modified frequency estimation (MFE) (Wang & Brennan, 2006) and 
       Braun-Holland under MFE (if requested).

   Input
     w1        weight for pop 1
     internal  internal anchor if non-zero; external anchor if zero
 	 nsv       number of score categories for v
     minv      minimum score for v
     maxv      maximum score for v      
     nsx       number of score categories for x
     minx      minimum score for x
     maxx      maximum score for x
     nsy       number of score categories for y
     miny      minimum score for y
     maxy      maximum score for y
     inc       increment--assumed to be the same for x, y, and v 
     bxvin[][] biv rel freq dist for x (rows) and v -- input 
     byvin[][] biv rel freq dist for y (rows) and v -- input 
     rv1       reliability for common items in pop 1 (used for MFE)
     rv2       reliability for common items in pop 2 (used for MFE)

   Output
     fxs[]      rel freq dist for x in syn pop  
     gys[]      rel freq dist for y in syn pop  
     eraw[]     equipercentile equated raw scores 
     a[]        slope for Braun-Holland Method 
     b[]        intercept for Braun-Holland Method 
     erawBH[]   Braun-Holland linear equated raw scores 

   NOTES: (1) Results are for FE if rv1 == 0 or rv2 == 0;
              otherwise results are for MFE
          (2) Assumes space allocated for fxs, gys, and eraw
              before function called
          (3) Braun-Holland results (a, b, erawBH) provided
              if a != NULL and b !=NULL and erawBH != NULL
          (4) If Braun-Holland results are desired, then before
              function is called, space must be allocated for erawBH;
              also, a and b variables must be declared and their
              addresses passed
          (5) Braun-Holland results are provided along with 
              equipercentile results; Braun Holland results not
              provided in isolation (refer to 2, above)

  Function calls other than C or NR utilities:
	 SyntheticDensities()
     cum_rel_freqs()
     perc_rank()
     EquiEquate()
     BH_LinEq()
                                                
  B. A. Hanson with updates by R. L. Brennan

  Date of last revision: 6/30/08 
*/
	int i,j;
	double **bxv,                    /* starts as bxvin; becomes matrix of x|v */
           **byv,                    /* starts as byvin; becomes matrix of y|v */
           *Fxs,                    /* cum relative freq dist for x in syn pop */
           *Gys,                    /* cum relative freq dist for y in syn pop */ 
           *prxs;                         /* percentile ranks for x in syn pop */

  bxv = dmatrix(0,nsx-1,0,nsv-1);
  byv = dmatrix(0,nsy-1,0,nsv-1);
  Fxs = dvector(0,nsx-1);
  Gys = dvector(0,nsy-1);
  prxs = dvector(0,nsx-1);

  /* initializations */

  for(j=0;j<nsx;j++) for(i=0;i<nsv;i++) bxv[j][i] = bxvin[j][i];
  for(j=0;j<nsy;j++) for(i=0;i<nsv;i++) byv[j][i] = byvin[j][i];
	
	/* fxs and gys: densities for x and y in synthetic population 
     Note: if 0 < rv1, rv2 <= 1, results are for MFE;
           otherwise, results are for FE */

  SyntheticDensities(w1,internal,nsv,nsx,bxv,nsy,byv,rv1,rv2,fxs,gys);

  /* Fxs, Gys, and prxs */

  cum_rel_freqs(minx,maxx,inc,fxs,Fxs); 
  cum_rel_freqs(miny,maxy,inc,gys,Gys);

  for (i=0;i<nsx;i++)
    prxs[i] = perc_rank(minx,maxx,inc,Fxs,score(i,minx,inc));
                         
  /* Equipercentile equating */     

  EquiEquate(nsy,miny,inc,Gys,nsx,prxs,eraw);

  /* Braun-Holland linear equating */

  if(a!=nullptr && b!=nullptr && erawBH!=nullptr){
    BH_LinEq(minx,maxx,miny,maxy,inc,fxs,gys,a,b);
    for(i=0;i<nsx;i++) erawBH[i] = *b + (*a)*(score(i,minx,inc));
  }

/* deallocations */                                                       

  free_dmatrix(bxv,0,nsx-1,0,nsv-1);
  free_dmatrix(byv,0,nsy-1,0,nsv-1);
  free_dvector(Fxs,0,nsx-1);
  free_dvector(Gys,0,nsy-1);
  free_dvector(prxs,0,nsx-1);
	
}

/******************************************************************************/

void SyntheticDensities(double w1, int internal, int nsv, int nsx, double **bxv, 
                        int nsy, double **byv, double rv1, double rv2,
                        double *fxs, double *gys)
/*
  Synthetic population densities for x and y

  Input
    w1       weight for population 1 (x)
    internal internal anchor if non-zero; external anchor if zero
    nsv      number of score categories for common items (v)
    nsx      number of score categories for x
    bxv      starts as bivariate rel freq dist for v (rows) and x (cols);
             ends as biv rel freq dist of x given v with v as rows
    nsy      number of score categories for y
    byv      starts as bivariate rel freq dist for y (rows) and v (cols);
             ends as biv rel freq dist of y given v with y as rows
    rv1      rel of common items in pop 1 (0-->FE; non-0-->MFE)
    rv2      rel of common items in pop 2 (0-->FE; non-0-->MFE)

  Output
    fxs    synthetic pop distribution (relative frequencies) for x
    gys    synthetic pop distribution (relative frequencies) for y

  NOTE: space for fxs and gys allocated before call

  Function calls other than C or NR utilities:
    MixSmooth()
    CondBivDist(() 
    ModCondBivDist()
    runerror()
                                                
  B. A. Hanson with updates by R. L. Brennan

  Date of last revision: 2/3/09 
*/
{
  double *h1,                /* rel freq dist for v in pop 1 after smoothing */
         *h2,                /* rel freq dist for v in pop 2 after smoothing */
         *f1,                 /* rel freq dist of x in pop 1 after smoothing */
         *g2,                 /* rel freq dist of y in pop 2 after smoothing */
         muv1=0.,                      /* mean of v in pop 1 after smoothing */
         muv2=0.;                      /* mean of v in pop 2 after smoothing */

	int i,v;
	
  h1 = dvector(0,nsv-1);                                      /* allocations */
  h2 = dvector(0,nsv-1);
  f1 = dvector(0,nsx-1);
  g2 = dvector(0,nsy-1);
			
  MixSmooth(nsv,nsx,1.0e-10,bxv,f1,h1); /* Smooth f1(x,v); get h1(v) & f1(x) */	
  MixSmooth(nsv,nsy,1.0e-10,byv,g2,h2); /* Smooth g2(y,v); get h2(v) & g2(y) */
 
	CondBivDist(nsx,nsv,bxv,h1);               /* get f1(x|v) with v as rows */
	CondBivDist(nsy,nsv,byv,h2);               /* get g2(y|v) with v as rows */
                                /* error corrected in avove statement 3-2-09 */

  /* modified conditional bivariate distributions for MFE */

  if(rv1 != 0 && rv2 != 0){

    if(rv1<0 || rv1>1 || rv2<0 || rv2>1) 
      runerror("\nrv1 or rv2 invalid");

    for(i=0;i<nsv;i++) muv1 += static_cast<double>(i)*h1[i];
    for(i=0;i<nsv;i++) muv2 += static_cast<double>(i)*h2[i];

    ModCondBivDist(internal,nsv,nsx,rv1,rv2,muv1,muv2,bxv);   /* mod f2(x|v) */
	ModCondBivDist(internal,nsv,nsy,rv2,rv1,muv2,muv1,byv);   /* mod g1(y|v) */

  }

  /* synthetic densities using Kolen and Brennan (2004, Eq 5.8) */  

	for(i=0;i<nsx;i++){
      fxs[i] = w1*f1[i];
      for(v=0;v<nsv;v++)
        fxs[i] += (1-w1)*bxv[i][v]*h2[v];              /* syn density for x */
    }

	for(i=0;i<nsy;i++){
      gys[i] = (1-w1)*g2[i];
      for(v=0;v<nsv;v++)
        gys[i] += w1*byv[i][v]*h1[v];                  /* syn density for y */
	}

  free_dvector(h1,0,nsv-1);                                /* deallocations */
  free_dvector(h2,0,nsv-1);
  free_dvector(f1,0,nsx-1);
  free_dvector(g2,0,nsy-1);

}

/****************************************************************************/

void MixSmooth(int nsv, int nsx, double unimix, double **bxv, 
               double *fx, double *hv)
/*
  "Smooth" bivariate distribution of v and x by mixing it with a 
  "small" uniform distribution.  The primary purpose of doing so is to 
  replace f(v)=0 with slightly positive relative frequencies.
  Note that descriptions and code are presented in terms of x and v,
  but function obviously applies to y and v, as well.
	
	Input
		nsv     number of categories for v
		nsx	    number of categories for x
		unimix	mixing proportion for uniform distribution; 1-unimix
				    is mixing proportion for observed distribution.
		bxv	    observed bivariate rel freqs for x (rows) and v (cols)

	Output
		bxv	    "smoothed" bivariate probabilities of x and v 
		fx		"smoothed" marginal univariate probabilities for x 
		hv		"smoothed" marginal univariate probabilities for v

    NOTE: space for fx and hv must be allocated before function called

  Function calls other than C or NR utilities: None
                                                
  B. A. Hanson with updates/revisions by R. L. Brennan

  Date of last revision: 6/30/08 
*/
{
  double srf,                /* rel freq "smoothed" so that it is never 0 */
	     obsmix,             /* weight for actual relative freq in a cell */
         uprob;                       /* cell probability = 1/(nrow*ncol) */
  int i,j;

  for(i=0;i<nsv;i++) hv[i] = 0.;
  for(j=0;j<nsx;j++) fx[j] = 0.;

  obsmix = 1.0 - unimix;     
  uprob = 1.0 / (static_cast<double>(nsv) * nsx);              
  unimix *= uprob;   /* unimix scaled so that sum of "smoothed" bxv[][]=1 */
  for(i=0;i<nsv;i++)
    for(j=0;j<nsx;j++){
      srf = (obsmix * bxv[j][i]) + unimix;
      bxv[j][i] = srf;                  /* "smoothed" elements of bxv[][] */
      hv[i] += srf;           /* "smoothed" rel freq of i-th element of v */
    }

  for(j=0;j<nsx;j++) 
    for(i=0;i<nsv;i++)    
      fx[j] += bxv[j][i];    /* "smoothed rel freq of j-th element of fx */
}

/**************************************************************************/

void CondBivDist(int nsx, int nsv, double **bxv, double *hv)
/*
	Conditional bivariate distribution of x given v. (Obviously applies
	as well to conditional bivariate distribution of y given v.)

	Input
		nsv    number of categoies for v
		nsx    number of categories for x 
 		bxv    bivariate relative freq distribution of x (rows) and v (cols)
		hv     mariginal rel freq distribution of v
 
	Output
		bxv	   conditional distribution of x given v (or y given v)
               where v is rows

  Function calls other than C or NR utilities: None
                                                
  B. A. Hanson with updates by R. L. Brennan

  Date of last revision: 6/30/08 
*/
{	
  int i,j;

  for(i=0;i<nsv;i++)
    for(j=0;j<nsx;j++)
      bxv[j][i] /= hv[i]; 
}

/************************************************************************/

void BH_LinEq(double minx, double maxx, double miny, double maxy, double inc, 
              double *fxs, double *gys, double *a, double *b)
/*
  Braun-Holland linear equating (based on FE or MFE systhetic densities

  Input
    minx  min score for x
    maxx  max score for x
    miny  min score for y
    maxy  max score for y
    inc   increment for both x and y
    fxs   synthetic density for x
    gsy   synthetic density for y

  Output
    a    slope
    b    intercept

  Function calls other than C or NR utilities:
    MomentsFromRFD()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08 
*/
{
  double mtsxs[4],                          /* moments for x in syn pop */
         mtsys[4];	                        /* moments for y in syn pop */
 
  MomentsFromRFD(minx,maxx,inc,nullptr,fxs,mtsxs); 
  MomentsFromRFD(miny,maxy,inc,nullptr,gys,mtsys);

  *a = mtsys[1]/mtsxs[1];                                      /* slope */
  *b = mtsys[0] - (*a)*mtsxs[0];                           /* intercept */

}

/**************************************************************************/

void ModCondBivDist(int internal, int nsv, int nsx, double rv1, double rv2, 
                    double muv1, double muv2, double **bxv)
/*
   For MFE method, modify conditional bivariate distribution of
   x given v with x as rows.  Description here worded 
   in terms of x given v, and code written in same manner, but 
   obviously function applies to y given v, as well.
 ... (comments omitted for brevity) ...
*/
{
  double v1, /* non integer v in pop 1 associated with integer v in pop 2 */
         slope = std::sqrt(rv2/rv1),                          /* slope for MFE */
         intercept = ((1-std::sqrt(rv2))/std::sqrt(rv1))*muv2 - 
                     ((1-std::sqrt(rv1))/std::sqrt(rv1))*muv1, /* intercept for MFE */
         **temp,                                 /* working values of bxv */
         **temp_collapsed,                       /* working values of buv */
         **temp_interpolated,               /* interpolated values of buv */
         sum;
  int v,j,u,
      v2,                                     /* integer v score in pop 2 */
      nsu;                       /* number of non-anchor score categories */
  
  temp = dmatrix(0,nsx-1,0,nsv-1);                            /* allocate */

  /****************************external anchor ****************************/
  if(internal==0){                                     
    for(v2=0;v2<nsv;v2++){                                 /* interpolate */
      v1 = intercept + slope*v2;             
      for(j=0;j<nsx;j++)    
        temp[j][v2] = interpolate(v1,nsv,bxv[j]);   
    }
  }
  /****************************internal anchor ****************************/
  else{                                          
    nsu = nsx - nsv + 1;
    temp_collapsed = dmatrix(0,nsu-1,0,nsv-1);                /* allocate */
    temp_interpolated = dmatrix(0,nsu-1,0,nsv-1);             /* allocate */

    /* store non-structural-zero elements of bxv[][] 
      in temp_collapsed[][]; in essense, collapse bxv[][] to buv[][] */

    for(v=0;v<nsv;v++){
      u = 0;
      for(j=0;j<nsx;j++){
        if(j<v || j>nsu-1+v) continue;           /* skip structural zeros */
        else temp_collapsed[u++][v] = bxv[j][v];
      }
    }

    /* get interpolated values under MFE for internal anchor */
  
    for(v2=0;v2<nsv;v2++){                              
      v1 = intercept + slope*v2;             
      for(j=0;j<nsu;j++)    
        temp_interpolated[j][v2] = interpolate(v1,nsv,temp_collapsed[j]);   
    }
  
    /* expand temp_interpolated[0...nsu][v2] to temp[0...nsx][v2] */
  
    for(v2=0;v2<nsv;v2++){
      u = 0;
      for(j=0;j<nsx;j++){               /* first, insert structural zeros */
        if(j<v2 || j>nsu-1+v2) temp[j][v2] = 0.; 
        else temp[j][v2] = temp_interpolated[u++][v2];
      }
    }

    free_dmatrix(temp_collapsed,0,nsu-1,0,nsv-1);           /* deallocate */
    free_dmatrix(temp_interpolated,0,nsu-1,0,nsv-1);        /* deallocate */
  }

  /* for both external and internal anchor
     normalize such that sum-over-x of f(x|v) = 1 */
    
  for(v=0;v<nsv;v++){     
    sum = 0.;
    for(j=0;j<nsx;j++) sum += temp[j][v];
    for(j=0;j<nsx;j++) 
      bxv[j][v] = (sum>0) ? temp[j][v] / sum : 0.;/* transfer to bxv[][] */
  }

  free_dmatrix(temp,0,nsx-1,0,nsv-1);  

}

/*************************************************************************/

void Chained_EE(int nsx, double *prx1, double minv, double maxv,
                double incv, int nsv, double *Gv1, double miny, 
                double incy, int nsy, double *Gy2, double *Gv2,
                double *eraw)
/*
  Chained equipercentile equating using the composed function
  discussed in Kolen and Brennan (2004, pp. 145-147)

  Input

    nsx    = number of score categories for X in pop 1
    prx1[] = PR for X in pop 1
    minv   = min score for V in both pop 1 and pop 2
    maxv   = max score for V in both pop 1 and pop 2
    incv   = increment for V in both pop 1 and pop 2
    nsv    = # of score categories for V in both pops
    Gv1[]  = crfd for V in pop 1

    miny   = min score for Y in pop 2
    incy   = increment for Y in pop 2
    nsy    = number of score categories for Y in pop 2
    Gy2[]  = crfd for Y in pop 2
    Gv2[]  = crfd for V in pop 2

  Output

    eraw[] = chained equipercentile Y-equivalents of X
    
  Notes.  (a) It is assumed that space for eraw[] allocated
              prior to function call
          (b) minv, maxv, incv, and nsv are 
              necessarily the same for both pops

  Function calls other than C or NR utilities:
    EquiEquate()
    perc_rank()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08 
*/
{
  int j;
  double *extov,                    /* V equivalents for X in pop 1 */
         *prv2;  /* PR's (wrt V for pop 2) for V equiv's in extov[] */
  
  extov = dvector(0,nsx-1);                             /* allocate */
  prv2 = dvector(0,nsx-1);   

  /* Put X on scale of V in pop 1; there are nsx
     (non-integer) V equivalents in extov[] */

  EquiEquate(nsv,minv,incv,Gv1,nsx,prx1,extov);

  /* Get PRs (relative to V for pop 2)
     for the non-integer V scores in extov[].
     Note that there are nsx equivalents */
  
  for (j=0;j<nsx;j++)   
    prv2[j] = perc_rank(minv,maxv,incv,Gv2,extov[j]);

  /* Using the PRs in prv2[] get the Y equivalents of X;
     i.e., put V equivalents on scale of Y */

  EquiEquate(nsy,miny,incy,Gy2,nsx,prv2,eraw);

  free_dvector(extov,0,nsx-1);                       /* deallocate */
  free_dvector(prv2,0,nsx-1);   
}

/*******************************************************************************/

void Print_SynDens(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
/*
  print synthetic densities under CINEG design along with FE or MFE results
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 
  
  Recall that methods are
    'M' = mean, 
    'L' = linear (Tucker, Levine-obs, Levine-true, chained)
    'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
    'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
    'G' = FE + BH-FE + MFE + BH-MFE
    'C' = Chained
    'A' = FE + BH-FE + MFE + BH-MFE + Chained

  Consequently, this function should be called only when 
    inall->method == 'E' || inall->method == 'F' || 
    inall->method == 'G' || inall->method == 'A'

  Function calls other than C or NR utilities: 
    MomentsFromRFD()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
  int i,j,k;
  double mtsfxs[4], mtsgys[4], fxssum, gyssum;

  if(inall->method == 'M' || inall->method == 'L' || inall->method == 'C')
    runerror("\nNot meaningful to call this function");  

  /* print syn densities and FE if inall->method == 'E', 'G' or 'A' */

  if(inall->method == 'E' || inall->method == 'G' || inall->method == 'A'){
    std::fprintf(fp,"\n\n%s\n\n",tt);
    std::fprintf(fp,"Synthetic Densities and Equated Scores Under\n"
               "Frequency Estimation (w1 = %7.5f)\n\n",inall->w1);
    std::fprintf(fp,"Input file for xv and pop 1: %s\n",inall->xfname);
    std::fprintf(fp,"Input file for yv and pop 2: %s\n\n",inall->yfname);
    std::fprintf(fp,"      Raw (x)           fxs");
    std::fprintf(fp,"      Raw (y)           gys");
    std::fprintf(fp,"          FE\n");

    for(i=0;i<inall->xv->ns1;i++){
      std::fprintf(fp,"\n %12.5f  ",score(i,inall->xv->min1,inall->xv->inc1));
      std::fprintf(fp,"%12.9f",r->fxs[0][i]);
      if(i<inall->yv->ns1){
        std::fprintf(fp," %12.5f  ",score(i,inall->yv->min1,inall->yv->inc1));
        std::fprintf(fp,"%12.9f",r->gys[0][i]);
      }
      else{
        std::fprintf(fp,"                            ");
      }
      std::fprintf(fp,"%12.5f",r->eraw[0][i]);
    }
    
    while(inall->yv->ns1>i){
      std::fprintf(fp,"               ");
      std::fprintf(fp,"\n %12.5f  ",score(i,inall->yv->min1,inall->yv->inc1));
      std::fprintf(fp,"%12.5f",r->gys[0][i]);
      i++;  
    }

    MomentsFromRFD(inall->xv->min1, inall->xv->max1, inall->xv->inc1,
                  nullptr, r->fxs[0], mtsfxs);
    MomentsFromRFD(inall->yv->min1, inall->yv->max1, inall->yv->inc1,
                  nullptr, r->gys[0], mtsgys);

    fxssum = gyssum = 0.;
    for (i=0;i<inall->xv->ns1;i++) fxssum += r->fxs[0][i];
    for (i=0;i<inall->yv->ns1;i++) gyssum += r->gys[0][i];

    std::fprintf(fp,"\n\n       Mean ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[0],mtsgys[0]);
    std::fprintf(fp,"\n       S.D. ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[1],mtsgys[1]);
    std::fprintf(fp,"\n       Skew ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[2],mtsgys[2]);;
    std::fprintf(fp,"\n       Kurt ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[3],mtsgys[3]);
    std::fprintf(fp,"\n\nSum of RF's ");
    std::fprintf(fp,"%15.9f            %15.9f",fxssum,gyssum);
  
    std::fprintf(fp,"\n\n");
    for(j=1;j<=63;j++) std::fprintf(fp,"*");
  }

  /* print syn densities and MFE if inall->method == 'F', 'G' or 'A' */

  if(inall->method == 'F' || inall->method == 'G' || inall->method == 'A'){
    std::fprintf(fp,"\n\n%s\n\n",tt);
    std::fprintf(fp,"Synthetic Densities and Equated Scores Under\n"
               "Modified Frequency Estimation (w1 = %7.5f)\n\n",inall->w1);
    std::fprintf(fp,"Input file for xv and pop 1: %s\n",inall->xfname);
    std::fprintf(fp,"Input file for yv and pop 2: %s\n\n",inall->yfname);
    std::fprintf(fp,"      Raw (x)           fxs");
    std::fprintf(fp,"      Raw (y)           gys");
    std::fprintf(fp,"         MFE\n");

    for(i=0;i<inall->xv->ns1;i++){
      std::fprintf(fp,"\n %12.5f  ",score(i,inall->xv->min1,inall->xv->inc1));
      std::fprintf(fp,"%12.9f",r->fxs[1][i]);
      if(i<inall->yv->ns1){
        std::fprintf(fp," %12.5f  ",score(i,inall->yv->min1,inall->yv->inc1));
        std::fprintf(fp,"%12.9f",r->gys[1][i]);
      }
      else{
        std::fprintf(fp,"                            ");
      }
      k = (inall->method == 'F') ? 0 : 2;
      std::fprintf(fp,"%12.5f",r->eraw[k][i]);
    }
    
    while(inall->yv->ns1>i){
      std::fprintf(fp,"               ");
      std::fprintf(fp,"\n %12.5f  ",score(i,inall->yv->min1,inall->yv->inc1));
      std::fprintf(fp,"%12.5f",r->gys[1][i]);
      i++;  
    }

    MomentsFromRFD(inall->xv->min1, inall->xv->max1, inall->xv->inc1,
                  nullptr, r->fxs[1], mtsfxs);
    MomentsFromRFD(inall->yv->min1, inall->yv->max1, inall->yv->inc1,
                  nullptr, r->gys[1], mtsgys);

    fxssum = gyssum = 0.;
    for (i=0;i<inall->xv->ns1;i++) fxssum += r->fxs[1][i];
    for (i=0;i<inall->yv->ns1;i++) gyssum += r->gys[1][i];

    std::fprintf(fp,"\n\n       Mean ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[0],mtsgys[0]);
    std::fprintf(fp,"\n       S.D. ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[1],mtsgys[1]);
    std::fprintf(fp,"\n       Skew ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[2],mtsgys[2]);;
    std::fprintf(fp,"\n       Kurt ");
    std::fprintf(fp,"%15.9f            %15.9f",mtsfxs[3],mtsgys[3]);
    std::fprintf(fp,"\n\nSum of RF's ");
    std::fprintf(fp,"%15.9f            %15.9f",fxssum,gyssum);
  
    std::fprintf(fp,"\n\n");
    for(j=1;j<=63;j++) std::fprintf(fp,"*");
  }
}
/*******************************************************************************/