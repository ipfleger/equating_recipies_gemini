/*
  Analytic_SEs.cpp 

  Analytic Standard Errors  

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
#include "NRutilities.hh"
#include "Analytic_SEs.h"
#include <cmath>

/*******************************************************************/

void SE_EPequate(int nsy, double npy, double *crfdy, int nsx, 
                 double incx, double npx, double *prdx, double *se_eeq)
/*
  Estimated standard errors for equipercentile equivalents 
  on scale of y (see Kolen & Brennan, 2004, p. 248). 
	
	Input
		nsy 	 Number of raw score categories for old form Y
		npy 	 Number of persons who took old form Y
		crfdy[]  Cumulative rel freq dist for old form Y
        nsx      Number of raw score categories for new form X
        incx     Increment between scores for x
        npx      Number of persons who took new form X
		prdx[]	 Percentile rank distribution for new form X 
 
	Output 
    se_eeq[]  standard error of equiperentile equivalents

    Note: assumes space for se_eeq[] allocated prior to call

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 3-17-09    


*/
{
  int i,j;
  double prp,                        /* pr expressed as proportion */
         g;                            /* relative frequency for y */

  for(j=0;j<=nsx-1;j++){
    prp = prdx[j]/100;
    if(prp>=1.) se_eeq[j] = 0.;
    else{
      for(i=1;i<=nsy-1;i++) if(crfdy[i] > prp) break;
      if(crfdy[i] != crfdy[i-1]){
        g = crfdy[i] - crfdy[i-1];
        se_eeq[j] = (1/(g*g))*
                    ( prp*(1-prp)*(npx+npy)/(npx*npy)
                      - (crfdy[i]-prp)*(prp-crfdy[i-1])/(npy*g) );
        se_eeq[j] =  (se_eeq[j] > 0) ? std::sqrt(se_eeq[j]) : 0.;

        /* The following line handles non-unit increments */
        se_eeq[j] = incx * se_eeq[j];
      } 
      else
        se_eeq[j] = 0.;
    }
  }

}

/*******************************************************************/