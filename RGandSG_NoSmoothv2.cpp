/*
  RGandSG_NoSmooth.cpp 

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


  File contains random groups and single group drivers and print function, along 
  with linear equating functions.  Equipercentile functions are in ERutilities.c                                             
*/

#include "ERutilities.h"
#include "NRutilities.h"
#include "RGandSG_NoSmooth.h"
#include <cstdio>
#include <cstring>


/*******************************************************************************/

void Wrapper_RN(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for getting mean, linear, or equipercentile equating for RG design
    with no smoothing

  Wrapper_RN() does the following:
      (i) allocate space for r->eraw[0][] and r->mts[0][] 
     (ii) populate elements of inall 
    (iii) compute r->a[0] = slope (if mean or linear), 
          r->b[0] = intercept (if mean or linear), 
          r->eraw[0][] and r->mts[0][] 
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep

  Input
  
    design = 'R' (random groups)
    method = 'M', 'L', or 'E'(mean, linear, or equipercentile)
    smoothing = 'N' (none)  
    x = struct USTATS
    y = struct USTATS 
    rep = replication number for bootstrap; should be set to 0
          for actual equating;  
    
  Output
    
    struct PDATA inall   populates selected values of inall 
    
    struct ERAW_RESULTS r          

      a[0] = slope (if method = 'M', 'L', or 'E')
      b[0] = intercept (if method = 'M', 'L', or 'E')
      eraw[][]: equated raw scores;          
                method (rows) by raw score (columns) matrix
                of equated scores. Here there is only one method.
                So, memory allocated for eraw[][] is: 
                eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                               (loc(x->max,x->min,x>-inc)]
                because we are getting equated raw scores for x   

      mts[][]:  moments for equated raw scores         
      
  NOTE: If Wrapper_RN() is called in a bootstrap loop, then in the 
        calling function struct ERAW_RESULTS must be different 
        from struct ERAW_RESULTS for the actual equating. 
                                            
  Function calls other than C or NR utilities:
    RGandSG_LinEq()
    EquiEquate()
    MomentsFromFD()  
                                                
  R. L. Brennan

  Date of last revision: 6/30/08 
*/
{ 
                          /* method names --- 10 characters; right justified */
  char *names[] ={"   y-equiv"};                    
  
  inall->rep = rep;                /* should be set to 0 for actual equating */
                     /* counting of replications done in Wrapper_Bootstrap() */ 
                    
  /* allocation and assignments for in
     Note that for every assignment of the form inall->(var) = xv->(var)
     values vary depending on whether x (or y) is for actual equating or
     a bootstrap sample; all other values are the same for the 
     actual equating and a bootstrap sample */
  
  if(inall->rep == 0){     /* no assignment or stor alloc for bootstrap reps */
    std::strcpy(inall->xfname,x->fname);
    std::strcpy(inall->yfname,y->fname);
    inall->x = x;
    inall->y = y;
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
    
    inall->nm = 1;                                                                    
    inall->names = cmatrix(0,inall->nm-1,0,11);    /* only one row/method, 0 */
    std::strcpy(inall->names[0],names[0]);
 
    inall->min = x->min;  
    inall->max = x->max;
    inall->inc = x->inc;
    inall->fdx = x->fd;
    inall->n = x->n;
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){    /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,0,0,loc(inall->max,inall->min,inall->inc)); 
    r->mts = dmatrix(0,0,0,3);            /* 0,3 is for the 4 moments */
  }
   
/* Compute equating results */
  
  if(method != 'E')
    RGandSG_LinEq(x->mts[0],x->mts[1],y->mts[0],y->mts[1],method,
                  inall->min,inall->max,inall->inc,r->a,r->b,r->eraw[0]);
  else
    EquiEquate(y->ns,y->min,y->inc,y->crfd,x->ns,x->prd,r->eraw[0]);

/* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

}

/*******************************************************************************/

void Wrapper_SN(char design, char method, char smoothing,  struct BSTATS *xy, 
               int rep, struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for getting mean, linear, or equipercentile equating for SG design
    with no smoothing

  Wrapper_SN() does the following:
      (i) allocate space for r->eraw[0][] and r->mts[0][] 
     (ii) populate elements of inall 
    (iii) compute r->a[0] = slope (if mean or linear),
          r->b[0] = intercept (if mean or linear), and
          r->eraw[0][] and r->mts[0][] 
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep
  NOTE: required input is one struct BSTATS

  Input
  
    design = 'S' (single group)
    method = 'M', 'L', or 'E'(mean, linear, or equipercentile)
    smoothing = 'N' (none)  
    xy = struct BSTATS 
    rep = replication number for bootstrap; should be set to 0
          for actual equating;  
    
    NOTE: it is assumed that the first variable
          in xy is indeed x and the second variable is y

  Output
    
    struct PDATA inall   populates selected values of inall 
    
    struct ERAW_RESULTS r          

      a[0] = slope (if mean or linear)
      b[0] = intercept (if mean or linear)
      eraw[][]: equated raw scores;          
                method (rows) by raw score (columns) matrix
                of equated scores. Here there is only one method.
                So, memory allocated for eraw[][] is: 
                eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                               (loc(x->max,x->min,x>-inc)]
                because we are getting equated raw scores for x   

      mts[][]:  moments for equated raw scores          
      
  NOTE: If Wrapper_SN() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 
                                            
  Function calls other than C or NR utilities:
    RGandSG_LinEq()
    EquiEquate()
    MomentsFromFD()  
                                                
  R. L. Brennan

  Date of last revision: 6/30/08   
*/
{ 
                       /* method names --- 10 characters; right justified */
  char *names[] ={"   y-equiv"};                    
  
  inall->rep = rep;                /* should be set to 0 for actual equating */
                   /* counting of replications done in Wrapper_Bootstrap() */ 
                    
  /* allocation and assignments for in
     Note that for every assignment of the form inall->(var) = xv->(var)
     values vary depending on whether xy is for actual equating or
     a bootstrap sample; all other values are the same for the 
     actual equating and a bootstrap sample */
  
  if(inall->rep == 0){     /* no assignment or stor alloc for bootstrap reps */
    std::strcpy(inall->xyfname,xy->fname);
    inall->xy = xy;
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
    
    inall->nm = 1;                                                                    
    inall->names = cmatrix(0,inall->nm-1,0,11);    /* only one row/method, 0 */
    std::strcpy(inall->names[0],names[0]);
 
    inall->min = xy->min1;  
    inall->max = xy->max1;
    inall->inc = xy->inc1;
    inall->fdx = xy->fd1;
    inall->n = xy->n;
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){         /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,0,0,loc(inall->max,inall->min,inall->inc)); 
    r->mts = dmatrix(0,0,0,3);                    /* 0,3 is for the 4 moments */
  }
   

/* Compute equating results */
  
  if(method != 'E')
    RGandSG_LinEq(xy->mts1[0],xy->mts1[1],xy->mts2[0],xy->mts2[1],method,
                  inall->min,inall->max,inall->inc,r->a,r->b,r->eraw[0]);
  else
    EquiEquate(xy->ns2,xy->min2,xy->inc2,xy->crfd2,xy->ns1,xy->prd1,r->eraw[0]);

/* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

}

/**********************************************************************/

void RGandSG_LinEq(double mnx, double sdx, double mny, double sdy,
                   char method, double min, double max, double inc, 
                   double *a, double *b, double *eraw)
/* Linear equating for RG and SG designs

  Input
     mnx = mean for x
     sdx = sd for x
     mny = mean for y
     sdy = sd for y
     method = 'M' or 'L'
     min = minimum score for x
     max = maximum score for x
     inc = increment for x

  Output
     a = slope
     b = intercept
     eraw[] = equated raw scores

  Function calls other than C or NR utilities:
    score() 
                                               
  R. L. Brennan

  Date of last revision: 6/30/08   
*/
{
  int j;

  *a = (method=='M') ? 1. : sdy/sdx;
  *b = mny - (*a)*mnx;  

/* get equated raw scores */               
  
   for(j=0;j<=loc(max,min,inc);j++)
      eraw[j] = *b + (*a)*score(j,min,inc); 

}

/*******************************************************************************/

void Print_RN(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
/*
  print RN results (RG; linear, mean or equipercentile; no smoothing)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
  int i,j;
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  if(inall->method == 'M') 
    std::fprintf(fp,"Mean Equating with Random Groups Design");
  else if(inall->method == 'L') 
    std::fprintf(fp,"Linear Equating with Random Groups Design");
  else  
    std::fprintf(fp,"Equipercentile Equating with Random Groups Design");
 
  std::fprintf(fp,"\n\nInput file for x: %s\n",inall->xfname);
  std::fprintf(fp,"Input file for y: %s\n\n",inall->yfname);
  
  if(inall->method != 'E'){        
    std::fprintf(fp,"      inter b  %12.5f\n",r->b[0]);
    std::fprintf(fp,"      slope a  %12.5f\n",r->a[0]);
  }

  for(i=1;i<=27;i++) std::fprintf(fp,"-");

 /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (x)  ");
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

}

/*******************************************************************************/

void Print_SN(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
/*
  print SN results (SG; linear, mean, or equipercentile ; no smoothing)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
  int i,j;
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  std::fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  if(inall->method == 'M') 
    std::fprintf(fp,"Mean Equating with Single Groups Design");
  else if(inall->method == 'L') 
    std::fprintf(fp,"Linear Equating with Single Groups Design");
  else  
    std::fprintf(fp,"Equipercentile Equating with Single Groups Design");
 
  std::fprintf(fp,"\n\nInput file for x and y: %s\n\n",inall->xyfname);
  
  if(inall->method != 'E'){        
    std::fprintf(fp,"      inter b  %12.5f\n",r->b[0]);
    std::fprintf(fp,"      slope a  %12.5f\n",r->a[0]);
  }
        
  for(i=1;i<=27;i++) std::fprintf(fp,"-");
  
 /* following code set up for any number of methods */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores");
  std::fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) std::fprintf(fp,"   Method %d:",j);
  std::fprintf(fp,"\nRaw Score (x)  ");
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
}

/*******************************************************************************/