/* IRT True Score Equating for Mixed-Format Tests

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

  For details about partial derivatives, refer to the following report:

  Kim, S., & Kolen, M.J. (2005). Methods for obtaining a common scale under
      unidimensional IRT models: A technical review and further extensions
      (Iowa Testing Programs Occasional Paper, No. 52). The University of
      Iowa.

  Note: All pointers are of 0-offset unless otherwise indicated.

*/

#include "IRTeq.h"
#include "IRTst.h"
#include "NRutilities.h" 
#include "ERutilities.h"

#include <cmath> // For std::pow

static struct IRTstControl *ContHandle; /* global setting to control method */
static struct ItemSpec *ContNewItems;
static double TrueS;

/***********************************************************************************/ 
short trueScoreEq(struct IRTstControl *Handle, struct ItemSpec *NewItems,
		struct ItemSpec *OldItems, int nScores, const double *newScores, 
		double *eqvOld, double *theta, double *newMin, double *OldMin) 
/*---------------------------------------------------------------------------*/		
{	
	int j, r, CatMax;
	long chance;
	double trueMinNew, trueMinOld, slope, intercept, xh, xl, x1, x2;
	double newScoreMin;	/* minimum new form score */
	double oldScoreMax;	/* Maximum old form score */
	double newTestMin;	/* Minimum new form test score */
	double oldTestMin;	/* Minimum old form test score */

	ContHandle   = Handle;
	ContNewItems = NewItems;

	trueMinNew = 0.0;
	trueMinOld = 0.0;
	oldScoreMax = 0.0;
	newTestMin = 0;
	oldTestMin = 0;

	for(j=1; j <= Handle->NewItemNum; j++) {
		if (NewItems[j].model == l3) {
			trueMinNew += (NewItems[j].ScoreFunc[1]  * (1.0 - NewItems[j].c[2])+ NewItems[j].ScoreFunc[2] * NewItems[j].c[2]);  /* added by Tianyou 7/18/08
/* the following original line by Seonghoon was commented out */
/*			trueMinNew += (NewItems[j].ScoreFunc[1] + (NewItems[j].ScoreFunc[2]) * (NewItems[j].c[2])); */
		} else {
			trueMinNew += NewItems[j].ScoreFunc[1];
		}

		newTestMin += NewItems[j].ScoreFunc[1];
	}

	for(j=1; j <= Handle->OldItemNum; j++) {
		if (OldItems[j].model == l3) {
			trueMinOld += (OldItems[j].ScoreFunc[1] + (OldItems[j].ScoreFunc[2]) * (OldItems[j].c[2]));
		} else {
			trueMinOld += OldItems[j].ScoreFunc[1];
		}
		CatMax = OldItems[j].CatNum;
		oldScoreMax += OldItems[j].ScoreFunc[CatMax];
		oldTestMin  += OldItems[j].ScoreFunc[1];
	}
	
	newScoreMin = *newScores;

	/* Find index of largest score below sum of c's */
	for (chance = 0; chance < nScores-1; chance++) {
		if (newScores[chance+1] > trueMinNew) break;
	}
	
	/* Compute slope and intercept of interpolating function for 
	   new form scores below chance */
	slope = trueMinNew - newTestMin;
	if (slope) slope = (trueMinOld - oldTestMin)/slope;
	intercept = oldTestMin - slope * newTestMin;

	xh=99; xl=-99;
	TrueS = newScores[chance+1];     /* Check if true score at theta=-99 is greater      */
	if (f_mix(xl) < 0.0) ++chance;   /* than computed largest score less than the        */
	                                 /* the sum of the c's, if so increase chance index. */

	/* Assign score equivalents */
	for(j=0; j <= chance; j++) {
		eqvOld[j] = intercept + newScores[j] * slope;
		theta[j] = xl;
	}

	x1 = -99.0;
	x2 = 99.0;

	for(j=chance+1; j < nScores-1; j++) {
		TrueS=newScores[j];
		for(r=0; r <= 10; r++) {
			theta[j] = er_rtsafe(funcd, x1, x2, 0.00001);
			if (theta[j] == -9999.0) theta[j] = -99.0;  /* tianyou added this line 7/18/08 */
			if(theta[j] != -9999.0 && theta[j] != -99.0 && theta[j] != 99.0) break;
			else {
				x1 = -std::pow(2.0, static_cast<double>(r));
				x2 =  std::pow(2.0, static_cast<double>(r));
			}
		}
		eqvOld[j]=trueScore(OldItems, Handle->OldItemNum, theta[j]);
	}
	
	/* Convert maximum score on new form to maximum score on old form */
	eqvOld[nScores-1] = oldScoreMax;
	theta[nScores-1] = xh;

	*newMin = trueMinNew;
	*OldMin = trueMinOld;
	
	return (0);		                                                  /* successful return */
}

/*****************************************************************************/

double f_mix(double theta)
/*
  Author: Seonghoon Kim
  Date of last revision 9/25/08
*/
{
	int j, k;
	double p, Wjk, v;
	
	v = 0.0;
	for(j=1; j <= ContHandle->NewItemNum; j++) {
		for(k=1; k <= ContNewItems[j].CatNum; k++) {
			p = ProbCCC(&ContNewItems[j], k, theta);
			Wjk = ContNewItems[j].ScoreFunc[k];
			v += Wjk*p;
		}
	}
	return (TrueS-v);
}

/*****************************************************************************/

double f_mixDer(double theta)
/*
  Author: Seonghoon Kim
  Date of last revision 9/25/08
*/
{
	int j, k;
	double pd, Wjk, v;
	
	v = 0.0;
	for(j=1; j <= ContHandle->NewItemNum; j++) {
		for(k=1; k <= ContNewItems[j].CatNum; k++) {
			pd = PdCCCoverTheta(&ContNewItems[j], k, theta);
			Wjk = ContNewItems[j].ScoreFunc[k];
			v += Wjk*pd;
		}
	}

	return -v;
}

/*****************************************************************************/
	
void funcd(double x, double *f, double *fd)
/*
  Author: Seonghoon Kim
  Date of last revision 9/25/08
*/
{
	*f = f_mix(x);
	*fd = f_mixDer(x);
	
	return;
}

/*****************************************************************************/

double trueScore(struct ItemSpec *Items, int n, double theta)
/*
  Author: Seonghoon Kim
  Date of last revision 9/25/08
*/
{
	int j, k;
	double p, Wjk, s;

	s = 0.0;
	for(j=1; j <= n; j++) {
		for(k=1; k <= Items[j].CatNum; k++) {
			p = ProbCCC(&Items[j], k, theta);
			Wjk = Items[j].ScoreFunc[k];
			s += Wjk*p;
		}
	}

	return (s);
}

/*****************************************************************************/

short IRTmixObsEq(struct IRTstControl *Handle, struct ItemSpec *NewItems,
		struct ItemSpec *OldItems, double wNew, double wOld,
		struct RawFitDist *newForm, struct RawFitDist *oldForm,
		struct RawTruObsEquiv *RawEq) 
/*------------------------------------------------------------------------------*/
{
	int i;
	int newMaxScrP, oldMaxScrP;
	double *newCumDist, *oldCumDist;
	double *newPR; /* percentile rank distribution for new form --Tianyou added on 8/18/08 */
	double oldFormMin, oldFormInc;
	
	newMaxScrP = newForm->nRaws;
	oldMaxScrP = oldForm->nRaws;
		
	/* fitted distribution for the new form with the new group */
	IRTmixObsDist(NewItems, Handle->NewItemNum, newMaxScrP, Handle->NewThetaNum,
					Handle->NewThetaValues, Handle->NewThetaWeights, &(newForm->nRaws),
					newForm->rawScrs, newForm->newFits);
	if(newForm->nRaws != newMaxScrP) 
		runerror("\nPossibly wrong execution 1 in IRTmixObsEq()\n");

	/* fitted distribution for the new form with the old group */
	IRTmixObsDist(NewItems, Handle->NewItemNum, newMaxScrP, Handle->OldThetaNum,
					Handle->OldThetaValues, Handle->OldThetaWeights, &(newForm->nRaws),
					newForm->rawScrs, newForm->oldFits);
	if(newForm->nRaws != newMaxScrP) 
		runerror("\nPossibly wrong execution 2 in IRTmixObsEq()\n");

	
	/* fitted distribution for the new form with the synthetic group */
	for(i=0; i < newForm->nRaws; i++) {
		newForm->synFits[i] = wNew*(newForm->newFits[i]) + wOld*(newForm->oldFits[i]);
	}
	
	/* fitted distribution for the old form with the new group */
	IRTmixObsDist(OldItems, Handle->OldItemNum, oldMaxScrP, Handle->NewThetaNum,
					Handle->NewThetaValues, Handle->NewThetaWeights, &(oldForm->nRaws),
					oldForm->rawScrs, oldForm->newFits);
	if(oldForm->nRaws != oldMaxScrP)
		runerror("\nPossibly wrong execution 3 in IRTmixObsEq()\n");

	/* fitted distribution for the old form with the old group */
	IRTmixObsDist(OldItems, Handle->OldItemNum, oldMaxScrP, Handle->OldThetaNum,
					Handle->OldThetaValues, Handle->OldThetaWeights, &(oldForm->nRaws),
					oldForm->rawScrs, oldForm->oldFits);
	if(oldForm->nRaws != oldMaxScrP)
		runerror("\nPossibly wrong execution 4 in IRTmixObsEq()\n");

	
	/* fitted distribution for the old form with the synthetic group */
	for(i=0; i < oldForm->nRaws; i++) {
		oldForm->synFits[i] = wNew*(oldForm->newFits[i]) + wOld*(oldForm->oldFits[i]);
	}
	
	newCumDist = static_cast<double *>(std::malloc(newForm->nRaws * sizeof(double)));
	if (!newCumDist) 
		runerror("memory allocation failure 1 in IRTmixObsEq()\n");

	newCumDist[0] = newForm->synFits[0];
	for(i=1; i < newForm->nRaws; i++) {
		newCumDist[i] = newCumDist[i-1] + newForm->synFits[i];
	}
	
	oldCumDist = static_cast<double *>(std::malloc(oldForm->nRaws * sizeof(double)));
	if (!oldCumDist) 
		runerror("memory allocation failure 2 in IRTmixObsEq()\n");

	oldCumDist[0] = oldForm->synFits[0];
	for(i=1; i < oldForm->nRaws; i++) {
		oldCumDist[i] = oldCumDist[i-1] + oldForm->synFits[i];
	}

	/* conduct equipercentile equating with synthetic distributions */
	oldFormMin = oldForm->rawScrs[0];
	oldFormInc = oldForm->rawScrs[1] - oldFormMin;

	newPR = static_cast<double *>(std::malloc(newForm->nRaws * sizeof(double)));
	if (!newPR) 
		runerror("memory allocation failure 3 in IRTmixObsEq()\n");


	/* compute the percentile rank function --Tianyou added 8/18/08 */
	for (i=0;i<newForm->nRaws;i++) {
	  newPR[i] = perc_rank(0, newForm->nRaws-1, 1, newCumDist, static_cast<double>(i));      /* row prd */
	}

	/* calling ERutilities equipercentile equating function (TW 8/18/08) */
	EquiEquate(oldForm->nRaws, oldFormMin, oldFormInc, oldCumDist,
               newForm->nRaws, newPR, RawEq->unroundedEqObs); 
	
	std::free(newCumDist);
	std::free(oldCumDist);
	std::free(newPR); // BUG FIX: Deallocate newPR
	
	return (0);  /* successful return */
}

/*****************************************************************************/

short IRTmixObsDist(struct ItemSpec *Items, int n, int MaxScrP, int nq, 
       double xqpts[], double xqwts[], int *nscr, double xscr[], double xmarg[])
/*------------------------------------------------------------------------------
  Functionality:
    Calculates the marginal distribution of total score for IRT models. The IRT
    models include the 3PL, LGR, GPC, and NR models.

  Input:
    Items:   pointer to designate items on a test form (0-offset)
    n :      number of items on a test form
    nq:      number of quadrature points
    xqpts:   vector for quadrature points (0-offset)
    xqwts:   vector for quadrature weights (0-offset)

  Output:
    nscr:    number of score categories for observed score distribution
    xscr:    vector of scores associated with each score category
             (are consecutive integers); 0-offset
    xmarg:   vector of marginal probabilities associated with each score
             category; 0-offset 

  Author: Seonghoon Kim
  Date of last revision 9/25/08

------------------------------------------------------------------------------*/
{
	int i, j, k, MaxCat;
	double theta, xx, *xnew;

	/* finds the maximum number of categories across items */
	MaxCat = 2;
	for(j=1; j <= n; j++) {
		if(Items[j].CatNum > MaxCat) MaxCat = Items[j].CatNum;
	}
	
	xnew = static_cast<double *>(std::malloc(MaxScrP * sizeof(double)));
	if (!xnew) 
		runerror("memory allocation failure 1 in IRTmixObsDist()\n");

	for (i=0; i < MaxScrP; i++) xmarg[i] = 0.0;

	for (k=1; k <= nq; k++) {
		theta = xqpts[k];
		if(ObsDistGivenTheta(theta, Items, n, MaxCat, MaxScrP, nscr, xscr, xnew)) {
			std::free(xnew);
			return (1); /* unsuccessful return */			
		}
		for (i=0; i < *nscr; i++) {
			xmarg[i] = xmarg[i] + xqwts[k]*xnew[i];
		}
	}

	xx = 0.0;
	for (i=0; i < *nscr; i++) xx += xmarg[i];
	for (i=0; i < *nscr; i++) xmarg[i] = xmarg[i] /xx;

	std::free(xnew);
	return (0);  /* successful return */
}

/*****************************************************************************/

short ObsDistGivenTheta(double theta, struct ItemSpec *Items, int n,
		int MaxCat, int MaxScrP, int *nscr, double xscr[], double xnew[])
/*------------------------------------------------------------------------------
  Functionality:
    Calculates the conditional distribution of total score given theta
    for IRT models. The IRT models include the 3PL, LGR, GPC, and NR models.
    Uses Hanson's (1994) generalization of the Lord-Wingersky (1982) recursive
    algorithm.

  Input:
    theta:   examinee ability
    Items:   pointer to designate items on a test form (0-offset)
    n :      number of items on a test form
    MaxCat:  the maximum number of categories across items
    MaxScrP: maximum number of score points

  Output:
    nscr:    number of score categories for observed score distribution
    xscr:    vector of scores associated with each score category
             (are consecutive integers); 0-offset
    xnew:    vector of probabilities associated with each score category;
             0-offset

  Author: Seonghoon Kim
  Date of last revision 9/25/08

------------------------------------------------------------------------------*/
{
	int i, j, k, index;
	int mino, maxo, minn, maxn;
	double *xitem,   /* zero-offset, but not use xitem[0] */
	       *xold;    /* zero-offset */

	xitem = static_cast<double *>(std::malloc((MaxCat+1)*sizeof(double)));
	if (!xitem) 
		runerror("memory allocation failure 1 in ObsDistGivenTheta()\n");

	
	xold = static_cast<double *>(std::malloc(MaxScrP * sizeof(double)));
	if (!xold) 
		runerror("memory allocation failure 2 in ObsDistGivenTheta()\n");


	/* calculates probabilities for Item 1 */
	for(k = 1; k <= Items[1].CatNum; k++) {
		xitem[k] = ProbCCC(&Items[1], k, theta);
	}

	mino = static_cast<int>(Items[1].ScoreFunc[1]);
	index = Items[1].CatNum;
	maxo = static_cast<int>(Items[1].ScoreFunc[index]);
	minn = mino;
	maxn = maxo;

	for(i=0; i <= maxn - minn; i++) {
		xold[i] = 0.0;
	}
	
	for(k=1; k <= Items[1].CatNum; k++) {
		index = static_cast<int>(Items[1].ScoreFunc[k]) - minn;
		xold[index] = xitem[k];  /* mino associated with index of 0 */
	}                           /* mino does vary; see below      */
	
	for(i=0; i <= maxn - minn; i++) {
		xnew[i] = xold[i];
	}

	if (n == 1) {
		for (i=0; i <= maxn - minn; i++) {
			xscr[i] = i + minn;
		}
		*nscr = maxn - minn + 1;
		std::free(xitem);
		std::free(xold);
		return (0); /* successful return */
	}

	/* updates distribution for items 2 through nitems */

	for (j=2; j <= n; j++) {
		for(k = 1; k <= Items[j].CatNum; k++) {
			xitem[k] = ProbCCC(&Items[j], k, theta);
		}
		
		if(recurs(mino, maxo, xold, Items[j].CatNum, Items[j].ScoreFunc, xitem, &minn, &maxn, xnew)) {
			std::free(xitem);
			std::free(xold);
			return (1);  /* unsuccessful return */
		}
		mino = minn;
		maxo = maxn;

		for (i=0; i <= maxn - minn; i++) {
			xold[i] = xnew[i];
		}
	}

	for (i=0; i <= maxn - minn; i++) {
		xscr[i] = i + minn;
	}
	*nscr = maxn - minn + 1;

	std::free(xitem);
	std::free(xold);
	return (0);  /* successful return */
}

/*****************************************************************************/

short recurs(int mino, int maxo, double xold[],
            int mitem, double iitem[], double xitem[],
            int * minn, int * maxn, double xnew[])
/*------------------------------------------------------------------------------
  Functionality:
    Updates a distribution of scores using Hanson's (1994) generalization of
    the Lord-Wingersky (1982) formula.

    Assumes that test scores are consecutive integers
    Assumes that item scores are consecutive integers
    Assumes that test and item scores are sorted from low to high.

  Input (r-1 added items):
    mino : minimum integer score for old distribution
    maxo : maximum integer score for old distribution
    xold : probability array for each score point from mino to maxo
          (f_r-1 (x|theta)); 0-offset
    mitem: number of distinct score points for the added (new) item
    iitem: values of the scoring function for the added item; 0-offset
    xitem: values of category response functions for the added item; 0-offset

  Output (r added items):
    minn : minimum integer score for new (updated) distribution
    maxn : maximum integer score for new (updated) distribution
    xnew : probability array for each score point from minn to maxn
          (f_r(x|theta)); 0-offset

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, j, in, io;

	*minn = mino + static_cast<int>(iitem[1]);
	*maxn = maxo + static_cast<int>(iitem[mitem]);

	for (i= *minn; i <= *maxn; i++) {
		in = i - *minn;
		xnew[in] = 0.0;
		for (j=1; j <= mitem; j++) {
			io = i - static_cast<int>(iitem[j]) - mino;
			if ( io >= 0 && io <= maxo-mino) {
				xnew[in] += xold[io]*xitem[j];
			}
		}
	}

	return (0);  /* successful return */
}

/*****************************************************************************/

double ProbCCC(struct ItemSpec *Item, int CatID, double theta)
/*------------------------------------------------------------------------------
  Functionality:
    Calculate the value of the category characteristic curve, P_jk, by model.

  Input:
    Item: One item of struct ItemSpec type 
    CatID: Response category ID in question (1 through CatNum)
    theta: ability value

  Output:
    Return the probability of an examinee having the value of theta
    responding in category ID
  Note:
    Uses arbitrarily the "old" scale with S = 1.0 and I = 0.0, since
    the functions in ProbDeriv.c are used.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double prob = 0.0;

	switch(Item->model)
	{
		case l3:
			prob = Prob3PL(CatID, theta, Item->ScaleConst,
				Item->a[2], Item->b[2], Item->c[2], "old", 1, 0);
			break;
		case gr:
			prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->a[2], Item->b, "old", 1, 0);
			break;
		case pc:
			prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->a[2], Item->b, "old", 1, 0);
			break;
		case nr:
			prob = ProbNRM(Item->CatNum, CatID, theta,
				Item->a, Item->c, "old", 1, 0);
			break;
		default:
			break;
	}
	return prob;
}

/*****************************************************************************/

double PdCCCoverTheta(struct ItemSpec *Item, int CatID, double theta)
/*------------------------------------------------------------------------------
  Functionality:
    By model, calculate the first (partial) derivative of P_jk
    with respect to ability theta.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pd = 0.0;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLoverTheta(CatID, theta, Item->ScaleConst,
				Item->a[2], Item->b[2], Item->c[2]);
			break;
		case gr:
			pd = PdLGRoverTheta(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->a[2], Item->b);
			break;
		case pc:
			pd = PdGPCoverTheta(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->a[2], Item->b);
			break;
		case nr:
			pd = PdNRMoverTheta(Item->CatNum, CatID, theta,
				Item->a, Item->c);
			break;
		default:
			break;
	}
	return pd;
}

/*****************************************************************************/

double Pd3PLoverTheta(int CatID, double theta, double D, double a, double b, double c)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate the first (partial) derivative of P_j
    with respect to theta.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double uj;

	uj = std::exp(D*a*(theta-b));
	if (CatID == 1) return ( -D*a*(1.0-c)*uj/((1.0+uj)*(1.0+uj)) );
	else return ( D*a*(1.0-c)*uj/((1.0+uj)*(1.0+uj)) );
}

/*****************************************************************************/

double PdLGRoverTheta(int CatNum, int CatID, double theta, double D, double a, double b[])
/*------------------------------------------------------------------------------
  Functionality:
    Under the GR model, calculate the first (partial) derivative of P_jk
    with respect to theta.
    
  Note:
    b[2] for the first item-step parameter
    b[CatNum] for the last item-step parameter

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double cp_jk, cp_jk1;

	if (CatID == 1) {
		cp_jk  = 1.0;
		cp_jk1 = 1.0/(1.0 + std::exp(-D*a*(theta-b[CatID+1])));
	}
	else {
		if (CatID < CatNum) {
			cp_jk  = 1.0/(1.0 + std::exp(-D*a*(theta-b[CatID])));
			cp_jk1 = 1.0/(1.0 + std::exp(-D*a*(theta-b[CatID+1])));
		}
		else { /* CatId == CatNum */
			cp_jk  = 1.0/(1.0 + std::exp(-D*a*(theta-b[CatID])));
			cp_jk1 = 0.0;
		}
	}
	return ( D*a*(cp_jk*(1.0-cp_jk) - cp_jk1*(1.0-cp_jk1)) );
}

/*****************************************************************************/

double PdGPCoverTheta(int CatNum, int CatID, double theta, double D, double a, double b[])
/*------------------------------------------------------------------------------
  Functionality:
    Under the GPC model, calculate the first (partial) derivative of P_jk
    with respect to ability theta.
    
  Note:
    b[2] for the actual first item-step parameter
    b[CatNum] for the last item-step parameter
    is dealt with as a special case of the nominal response model

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k, l;
	double ujk, vj = 0.0;
	double a_sum, b_sum, a_u_sum = 0.0;

	for (k = 1; k <= CatNum ; k++) {
		a_sum = static_cast<double>(k)*D*a;
		b_sum = 0.0;
		for (l = 2; l <= k ; l++) {  /* b[1] = 0 */
			b_sum += b[l];
		}
		b_sum *= -D*a;
		ujk = std::exp(a_sum*theta + b_sum);
		vj += ujk;
		a_u_sum += a_sum*ujk;
	}

	a_sum = static_cast<double>(CatID)*D*a;
	b_sum = 0.0;
	for (l = 2; l <= CatID ; l++) {
		b_sum += b[l];
	}
	b_sum *= -D*a;
	ujk = std::exp(a_sum*theta + b_sum);
	
	return ( a_sum*ujk/vj - ujk*a_u_sum/(vj*vj) );
}

/*****************************************************************************/

double PdNRMoverTheta(int CatNum, int CatID, double theta, double a[], double c[])
/*------------------------------------------------------------------------------
  Functionality:
    Under the NR model, calculate the first (partial) derivative of P_jk
    with respect to ability theta.
    
  Note:
    a[1..CatNum] for the discrimination parameters
    c[1..CatNum] for the intercept parameters
    No scaling constant

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double ujk, vj = 0.0;
	double a_u_sum = 0.0;

	for (k = 1; k <= CatNum ; k++) {
		ujk = std::exp(a[k]*theta + c[k]);
		vj += ujk;
		a_u_sum += a[k]*ujk;
	}
	
	ujk = std::exp(a[CatID]*theta + c[CatID]);

	return ( a[CatID]*ujk/vj - ujk*a_u_sum/(vj*vj) );
}

/*****************************************************************************/

void RawFitMem(struct ItemSpec *Items, const char *oldOrnew, 
		struct IRTstControl *Handle, struct RawFitDist *Form)
/*------------------------------------------------------------------------------*/
{
	int i, j, MinScrP,                              /* minimum possible score */
              MaxScrP,    /* max number of score categories --- not max score */
              ItemNum, index;
	double temp1=0, temp2=0;

	if (std::strcmp(oldOrnew, "old")==0) ItemNum = Handle->OldItemNum;
	else                            ItemNum = Handle->NewItemNum;

	/* Calculate the possible largest number of observed score categories */
	
	MaxScrP = MinScrP = 0;

	for(j=1; j <= ItemNum; j++) {
		index = Items[j].CatNum;
		temp1 += Items[j].ScoreFunc[1]; 
		temp2 += Items[j].ScoreFunc[index];
	}
	if (temp1<0) 
		MinScrP = static_cast<int>(temp1 - .5);
	else 
		MinScrP = static_cast<int>(temp1 + .5);  /* +.5: rounding the double to integer */

	MaxScrP = static_cast<int>(temp2) - MinScrP + 1;      /* maximum number of score categories */

	if (std::strcmp(oldOrnew, "old")==0) {
	  Handle->OldRawMin = MinScrP;
	  Handle->OldRawMax = temp2;
	  Handle->OldRawInc = 1;
	}
	else {
      Handle->NewRawMin = MinScrP;
	  Handle->NewRawMax = temp2;
	  Handle->NewRawInc = 1;
	}

	Form->nRaws = MaxScrP;               /* maximum number of score categories */
	
	Form->rawScrs = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(Form->rawScrs)) 
		runerror("memory allocation failure 1 in EqRawFitMem()\n");

	Form->newFits = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(Form->newFits)) 
		runerror("memory allocation failure 2 in EqRawFitMem()\n");

	Form->oldFits = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(Form->oldFits))
		runerror("memory allocation failure 3 in EqRawFitMem()\n");

	Form->synFits = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(Form->synFits)) 
		runerror("memory allocation failure 4 in EqRawFitMem()\n");
 	
	for(i=0; i < MaxScrP; i++) {
		Form->rawScrs[i] = MinScrP + i;
	}

	return;
}

/*****************************************************************************/

void RawEqResultsMem(struct ItemSpec *NewItems, struct IRTstControl *Handle,
		struct RawTruObsEquiv *RawEq)
/*------------------------------------------------------------------------------*/
{
	int j, MinScrP, MaxScrP, ItemNum, index;
	double temp1=0, temp2=0;

	ItemNum = Handle->NewItemNum;

	/* Calculate the possible largest number of observed score categories */

	MaxScrP = MinScrP = 0;

	for(j=1; j <= ItemNum; j++) {
		index = NewItems[j].CatNum;
		temp1 += NewItems[j].ScoreFunc[1]; 
		temp2 += NewItems[j].ScoreFunc[index];
	}
	if (temp1<0) 
		MinScrP = static_cast<int>(temp1 - .5);
	else 
		MinScrP = static_cast<int>(temp1 + .5);  /* +.5: rounding the double to integer */

	MaxScrP = static_cast<int>(temp2) - MinScrP + 1;
	
	RawEq->nRaws = MaxScrP;
	
	RawEq->thetaTru = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(RawEq->thetaTru)) 
		runerror("memory allocation failure 1 in EqResultMem()\n");

	RawEq->unroundedEqTru = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(RawEq->unroundedEqTru)) 
		runerror("memory allocation failure 2 in EqResultMem()\n");

	RawEq->unroundedEqObs = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(RawEq->unroundedEqObs)) 
		runerror("memory allocation failure 3 in EqResultMem()\n");

	RawEq->roundedEqTru = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(RawEq->roundedEqTru)) 
		runerror("memory allocation failure 4 in EqResultMem()\n");

	RawEq->roundedEqObs = static_cast<double *>(std::malloc(MaxScrP*sizeof(double)));
	if (!(RawEq->roundedEqObs)) 
		runerror("memory allocation failure 5 in EqResultMem()\n");

	return;
}

/*****************************************************************************/

void RawFitDeAlloc(struct RawFitDist *Form)
/*------------------------------------------------------------------------------*/
{
	std::free(Form->rawScrs);
	std::free(Form->newFits);
	std::free(Form->oldFits);
	std::free(Form->synFits);
	
	return;
}

/*****************************************************************************/

void RawEqResultsDeAlloc(struct RawTruObsEquiv *RawEq)
/*------------------------------------------------------------------------------*/
{
	std::free(RawEq->thetaTru);
	std::free(RawEq->unroundedEqTru);
	std::free(RawEq->unroundedEqObs);
	std::free(RawEq->roundedEqTru);
	std::free(RawEq->roundedEqObs);
	
	return;
}

/**************************************************************************************/

 void Wrapper_IRTeq(char design, char method, double w1, 
     char ItemNewFile[], char ItemOldFile[], char DistNewFile[], char DistOldFile[], 
     struct ItemSpec *NewItems, struct ItemSpec *OldItems, 
     struct RawFitDist *NewForm, struct RawFitDist *OldForm,
     struct RawTruObsEquiv *RawEq, struct IRTstControl *StInfo, int *NewFD, 
     struct IRT_INPUT *irtall, struct PDATA *pinall, struct ERAW_RESULTS *r)              
/*------------------------------------------------------------------------------*/
{ 
  FILE *inf;
  int i;                                   
  char *names[] ={"IRT Tr Scr", "IRT Ob Scr"};                /* method names */

  /* read the input for the new form items */

  inf = std::fopen(ItemNewFile, "r");
  if(!inf) runerror("can't open file for new form items \n");
  NewItems = ItemInfoRead(inf, "new", StInfo);
  std::fclose(inf);

  /* read the input for the old form items */

  inf = std::fopen(ItemOldFile, "r");
  if(!inf) runerror("can't open file for old form items \n");
  OldItems = ItemInfoRead(inf, "old", StInfo);
  std::fclose(inf);  

  /* The following code segment reads quadrature distributions
     only if IRT observed score equating is requested, since
     quadrature distributions are not required for IRT
     true score equating */

  if(method == 'O' || method == 'A'){

    /* read the input for the new group's ability distribution */

    inf = std::fopen(DistNewFile, "r");
    if(!inf) runerror("can't open file for new form ability distribution \n");
    ThetaInfoRead(inf, "new", StInfo);
    std::fclose(inf);

    /* read the input for the old group's ability distribution */

    inf = std::fopen(DistOldFile, "r");
    if(!inf) runerror("can't open file for old form ability distribution \n");
    ThetaInfoRead(inf, "old", StInfo);
    std::fclose(inf);
  }
  
  /* Allocate memory */

  RawEqResultsMem(NewItems, StInfo, RawEq);
  RawFitMem(OldItems, "old", StInfo, OldForm); 
  RawFitMem(NewItems, "new", StInfo, NewForm);

  /* assignments for irtall */
  
  irtall->method = method;
  irtall->stControl = *StInfo;
  irtall->NewItems = NewItems;
  irtall->OldItems = OldItems;
  irtall->RawEq = RawEq;
  irtall->NewForm = NewForm;
  irtall->OldForm = OldForm;

  std::strcpy(irtall->ItemNewFile,ItemNewFile);
  std::strcpy(irtall->ItemOldFile,ItemOldFile); 
  if(method == 'O' || method == 'A'){ 
    std::strcpy(irtall->DistNewFile,DistNewFile);
    std::strcpy(irtall->DistOldFile,DistOldFile);
  }

  /* assignments for pinall */

  pinall->rep = 0;                            /* required for Wrapper_ESS() */
  pinall->design = design;
  pinall->w1 = w1;
  pinall->min = StInfo->NewRawMin;
  pinall->max = StInfo->NewRawMax;
  pinall->inc = StInfo->NewRawInc;
  pinall->names = cmatrix(0,1,0,11);                           /* two names */
 
  if(method == 'T'){                         /* method = true-score equating */
    pinall->nm = 1;
    for(i=0;i<1;i++) std::strcpy(pinall->names[i],names[i]);
  }
  else if(method == 'O'){                 /* method = observed score equating*/
    pinall->nm = 1;
    for(i=1;i<2;i++) std::strcpy(pinall->names[i],names[i]);
  }
  else if(method == 'A'){             /* method = both true and obs equating */
    pinall->nm = 2;
    for(i=0;i<2;i++) std::strcpy(pinall->names[i],names[i]);
  }
  else runerror("\nInvalid method\n");

  /* allocation and assignments for r (equivalents and moments) */ 

  r->eraw = dmatrix(0,pinall->nm-1,0,NewForm->nRaws-1);                
  r->mts = dmatrix(0,pinall->nm-1,0,3);          /* 0,3 is for the 4 moments */                            
   
  /* IRT true score equating */
                                           
  if(method == 'T') {                     
    if(trueScoreEq(StInfo, NewItems, OldItems, NewForm->nRaws, NewForm->rawScrs,
           RawEq->unroundedEqTru, RawEq->thetaTru,
           &(RawEq->TCCnewMin), &(RawEq->TCColdMin))) 
      runerror("\nIRT true score equating failed!\n");
    for(i=0;i<NewForm->nRaws;i++) 
      r->eraw[0][i] = RawEq->unroundedEqTru[i];
  }

  /* IRT observed score equating.  First code segment declares an error and exits
     if IRT observed score equating cannot be done because 
     w1<0 || w1>1 || DistNewFile == NULL || DistOldFile == NULL*/

  if(method == 'O' || method =='A') 
    if(w1<0 || w1>1 || std::strcmp(DistNewFile,"NULL") == 0 || std::strcmp(DistOldFile,"NULL") == 0)
      runerror("\nIRT observed-score equating cannot be performed because\n"
               "w1<0 || w1>1 || DistNewFile == NULL || DistOldFile == NULL \n");

  if(method == 'O') {
    if(IRTmixObsEq(StInfo, NewItems, OldItems, w1, 1-w1, NewForm, OldForm, RawEq)) 
      runerror("\nIRT observed score equating failed!\n");
    for(i=0;i<NewForm->nRaws;i++) 
      r->eraw[0][i] = RawEq->unroundedEqObs[i];
  }

  /* Both IRT true-score and observed-score equating */

  if(method == 'A'){                /* IRT true + IRT observed score equating */
    if(trueScoreEq(StInfo, NewItems, OldItems, NewForm->nRaws, NewForm->rawScrs,
           RawEq->unroundedEqTru, RawEq->thetaTru,
           &(RawEq->TCCnewMin), &(RawEq->TCColdMin)))
      runerror("\nIRT true score equating failed!\n");

    if(IRTmixObsEq(StInfo, NewItems, OldItems, w1, 1-w1, NewForm, OldForm, RawEq))
      runerror("\nIRT observed score equating failed!\n"); 
    for(i=0;i<NewForm->nRaws;i++) {
      r->eraw[0][i] = RawEq->unroundedEqTru[i];
      r->eraw[1][i] = RawEq->unroundedEqObs[i];
    }
  }
                        
/* get moments based on actual frequencies for group that took new form*/

  if (NewFD!=nullptr) {
    irtall->NewFD = NewFD;
    pinall->fdx = NewFD;
    for(i=0;i<=pinall->nm-1;i++) 
      MomentsFromFD(StInfo->NewRawMin,StInfo->NewRawMax,StInfo->NewRawInc,
                    r->eraw[i],NewFD,r->mts[i]);
  }
  else {
    irtall->NewFD = nullptr;
    pinall->fdx = nullptr;
  }

  /* calculate raw score moments for both new form and old form
     using IRT fitted distributions for new, old, and synthetic gps
     and the quadrature distributions for the new and old groups.
     To get these results, the quadrature distributions must be
     provided, as well as w1.  These results have an ambiguous status 
     in the context of IRT true-score equating which does not 
     involve synthetic groups. 

     From here to the end of the function there are 8 calls to MomentsFromRFD().
     The first three parameters of each call were revised on 3/8/09  */

  if(method == 'O' || method == 'A'){

    MomentsFromRFD(StInfo->NewRawMin,StInfo->NewRawMax,StInfo->NewRawInc, 
                   NewForm->rawScrs, NewForm->newFits, NewForm->mtsng);

    MomentsFromRFD(StInfo->NewRawMin,StInfo->NewRawMax,StInfo->NewRawInc, 
                   NewForm->rawScrs, NewForm->oldFits, NewForm->mtsog);

    MomentsFromRFD(StInfo->NewRawMin,StInfo->NewRawMax,StInfo->NewRawInc,
                   NewForm->rawScrs, NewForm->synFits, NewForm->mtssg);

    MomentsFromRFD(StInfo->OldRawMin,StInfo->OldRawMax,StInfo->OldRawInc,
                   OldForm->rawScrs, OldForm->newFits, OldForm->mtsng);

    MomentsFromRFD(StInfo->OldRawMin,StInfo->OldRawMax,StInfo->OldRawInc,
                   OldForm->rawScrs, OldForm->oldFits, OldForm->mtsog);

    MomentsFromRFD(StInfo->OldRawMin,StInfo->OldRawMax,StInfo->OldRawInc,
                   OldForm->rawScrs, OldForm->synFits, OldForm->mtssg);
  }

  /* moments for true-score and observed-score equivalents 
     using quadrature distribution for group that took new form.
     These differ from the GroupFormMoments() in that these
     employ the equivalents. Caveat: the results are somewhat 
     ambiguous for IRT true-score equating which does not i
     involve synthetic groups */

  if(method == 'T' || method == 'A')
    MomentsFromRFD(StInfo->NewRawMin,StInfo->NewRawMax,StInfo->NewRawInc, 
                   RawEq->unroundedEqTru, NewForm->newFits, RawEq->mtsTru);

  if(method == 'O' || method == 'A')
    MomentsFromRFD(StInfo->NewRawMin,StInfo->NewRawMax,StInfo->NewRawInc, 
                   RawEq->unroundedEqObs, NewForm->newFits, RawEq->mtsObs);

  pinall->IRT_Input = irtall;
 
  return;
} 

/*********************************************************************************/

void Print_IRTeq(FILE *fp, char tt[], struct PDATA *pinall, 
                 struct ERAW_RESULTS *r, int PrintFiles)
/*
  print Wrapper_IRTeq results (IRT true score, observed score equating, or both)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    pinall =  struct PDATA
    r = struct ERAW_RESULTS 
    PrintFiles: 0 --> don't print item parameter and quadrature files
                1 --> print item parameter and quadrature files

  Function calls other than C or NR utilities: None
                                                
  Authors: T. D. Wang and R. L. Brennan
  Date of last revision 9/15/08
*/
{
  int i,j, 
      xmax = loc(pinall->IRT_Input->stControl->NewRawMax,
                 pinall->IRT_Input->stControl->NewRawMin,
                 pinall->IRT_Input->stControl->NewRawInc); /* loc of max X score */
  double moments[4];

  /* Print introductory information */
  
  std::fprintf(fp,"\n\n%s\n\n",tt);
    
  if (pinall->IRT_Input->method=='T')
    std::fprintf(fp,"IRT True Score Equating\n\n");
  else if (pinall->IRT_Input->method=='O')
    std::fprintf(fp,"IRT Observed Score Equating\n\n");
  else if (pinall->IRT_Input->method=='A')
    std::fprintf(fp,"IRT True Score Equating and Observed Score Equating");

  std::fprintf(fp,"\n\nInput Design = %c (ignored)",pinall->design);
  std::fprintf(fp,"\n\nName of file containing item parameter estimates for X: %s",
    pinall->IRT_Input->ItemNewFile);  
  if(PrintFiles == 1) Print_file(pinall->IRT_Input->ItemNewFile,fp);

  std::fprintf(fp,"\nName of file containing item parameter estimates for Y: %s",
    pinall->IRT_Input->ItemOldFile);
  if(PrintFiles == 1) Print_file(pinall->IRT_Input->ItemOldFile,fp);

  if(pinall->IRT_Input->method=='O' || pinall->IRT_Input->method=='A'){
    std::fprintf(fp,"\n\nw1 = weight for new group that took X = %7.5f\n",pinall->w1);

    std::fprintf(fp,"\nName of file containing theta distribtion for X: %s",
      pinall->IRT_Input->DistNewFile);
    if(PrintFiles == 1) Print_file(pinall->IRT_Input->DistNewFile,fp);

    std::fprintf(fp,"\nName of file containing theta distribtion for Y: %s",
      pinall->IRT_Input->DistOldFile);
    if(PrintFiles == 1) Print_file(pinall->IRT_Input->DistOldFile,fp);
  }

  std::fprintf(fp,"\n\nNumber of score categories for the new form  = %d\n",
           pinall->IRT_Input->NewForm->nRaws);
  std::fprintf(fp,"Number of score categories for the old form  = %d\n",
           pinall->IRT_Input->OldForm->nRaws);

  std::fprintf(fp,"\n");
  for(i=1;i<=56;i++) std::fprintf(fp,"-");

 /* Print equivalents (code set up for any number of methods).
    Also, print theta used to obtain IRT true score equivalents,
    if IRT true score equating is performed */

  std::fprintf(fp,"\n\n");
  for(j=1;j<=17+(pinall->nm*12-18)/2;j++) std::fprintf(fp," ");
  std::fprintf(fp,"Equated Raw Scores" );
  std::fprintf(fp,"\n\n              ");
  for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"    Method %d:",j);
  if(pinall->IRT_Input->method == 'T' || pinall->IRT_Input->method == 'A') 
    std::fprintf(fp,"       theta for");
  std::fprintf(fp,"\nRaw Score (X)  ");

  if(pinall->IRT_Input->method == 'O') std::fprintf(fp,"  %s",pinall->names[1]);
  else
    for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"  %s",pinall->names[j]);

  if(pinall->IRT_Input->method == 'T' ||pinall->IRT_Input->method == 'A') 
    std::fprintf(fp,"       %s",pinall->names[0]);
  
  std::fprintf(fp,"\n\n");

  for(i=0;i<=xmax;i++){
    std::fprintf(fp,"\n %12.5f  ",score(i,pinall->min,pinall->inc));
    for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"%12.5f",r->eraw[j][i]);
    if(pinall->IRT_Input->method == 'T' || pinall->IRT_Input->method == 'A')
      std::fprintf(fp,"     %12.5f",pinall->IRT_Input->RawEq->thetaTru[i]); 
  }

  /* print moments based on actual FD for group 1 that took Form X */

  if (pinall->IRT_Input->NewFD!=nullptr) {
    std::fprintf(fp,"\n\n");
    for(i=1;i<=56;i++) std::fprintf(fp,"-");
    std::fprintf(fp,"\n\n");
    std::fprintf(fp, "Moments based on actual frequency dist for new form\n");
    std::fprintf(fp,"\n\n         Mean  ");
    for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][0]);
    std::fprintf(fp,"\n         S.D.  ");
    for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][1]);
    std::fprintf(fp,"\n         Skew  ");
    for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][2]);
    std::fprintf(fp,"\n         Kurt  ");
    for(j=0;j<=pinall->nm-1;j++) std::fprintf(fp,"%12.5f",r->mts[j][3]);
  }

  /* If IRT observed score equating is requested, print moments based on 
     quadrature distribution for group that took Form X.  Note that
     if method == 'A', the moments for IRT true score equating are
     reported, even though IRT true score equating is group invariant,
     at least in theory */

  if(pinall->IRT_Input->method == 'O' || pinall->IRT_Input->method == 'A'){ 
    std::fprintf(fp,"\n\n");
    for(i=1;i<=56;i++) std::fprintf(fp,"-");
    std::fprintf(fp,"\n\n");	              /* change made in next line 3-08-09 */
    std::fprintf(fp, "Moments based on IRT fitted dist for new form\n");  
    std::fprintf(fp,"\n\n         Mean  ");
    if(pinall->IRT_Input->method == 'O') 
      std::fprintf(fp,"%12.5f",pinall->IRT_Input->RawEq->mtsObs[0]);
    else std::fprintf(fp,"%12.5f%12.5f",pinall->IRT_Input->RawEq->mtsTru[0],
                                   pinall->IRT_Input->RawEq->mtsObs[0]);
    std::fprintf(fp,"\n         S.D.  ");
    if(pinall->IRT_Input->method == 'O') 
      std::fprintf(fp,"%12.5f",pinall->IRT_Input->RawEq->mtsObs[1]);
    else std::fprintf(fp,"%12.5f%12.5f",pinall->IRT_Input->RawEq->mtsTru[1],
                                   pinall->IRT_Input->RawEq->mtsObs[1]);
    std::fprintf(fp,"\n         Skew  ");
    if(pinall->IRT_Input->method == 'O') 
      std::fprintf(fp,"%12.5f",pinall->IRT_Input->RawEq->mtsObs[2]);
    else std::fprintf(fp,"%12.5f%12.5f",pinall->IRT_Input->RawEq->mtsTru[2],
                                   pinall->IRT_Input->RawEq->mtsObs[2]);
    std::fprintf(fp,"\n         Kurt  ");
    if(pinall->IRT_Input->method == 'O') 
      std::fprintf(fp,"%12.5f",pinall->IRT_Input->RawEq->mtsObs[3]);
    else std::fprintf(fp,"%12.5f%12.5f",pinall->IRT_Input->RawEq->mtsTru[3],
                                   pinall->IRT_Input->RawEq->mtsObs[3]);
  }

  std::fprintf(fp,"\n\n");
  for(i=1;i<=56;i++) std::fprintf(fp,"-");
  std::fprintf(fp,"\n\n");

  /* If IRT observed score equating is performed, print distributions
     for old and new forms for group 1, group 2, and synthetic group.
     Also print moments. All of these distributions are smoothed in the 
     IRT sense, i.e., they are expected number-of-points distributions given
     the quadrature distributions.  This means, for example, that the 
     group 1 moments for X reported here will not generally be the 
     same as the moments based on the actual FD for group 1.  */

  if(pinall->IRT_Input->method == 'O' || pinall->IRT_Input->method == 'A'){

    std::fprintf(fp, "Relative Frequency Distributions for New Form X\n"
		        "          Using IRT Fitted Distributions\n\n");

    std::fprintf(fp, "Raw Score (X)       grp 1       grp 2       grp s\n");
    for(i=0; i < pinall->IRT_Input->NewForm->nRaws; i++) 
      std::fprintf(fp, "   %10.5f  %10.5f  %10.5f  %10.5f\n",
        pinall->IRT_Input->NewForm->rawScrs[i], pinall->IRT_Input->NewForm->newFits[i],
        pinall->IRT_Input->NewForm->oldFits[i], pinall->IRT_Input->NewForm->synFits[i]);

    std::fprintf(fp, "\n\nRelative Frequency Distributions for Old Form Y\n"                   
		            "          Using IRT Fitted Distributions\n\n");

    std::fprintf(fp, "Raw Score (Y)       grp 1       grp 2       grp s\n");
    for(i=0; i < pinall->IRT_Input->OldForm->nRaws; i++)
      std::fprintf(fp, "   %10.5f  %10.5f  %10.5f  %10.5f\n",
        pinall->IRT_Input->OldForm->rawScrs[i], pinall->IRT_Input->OldForm->newFits[i],
        pinall->IRT_Input->OldForm->oldFits[i], pinall->IRT_Input->OldForm->synFits[i]);

    std::fprintf(fp,"\n");
    for(i=1;i<=56;i++) std::fprintf(fp,"-");

    /* Raw score moments (based on quadrature distributions)*/
                                       /* change in next line made on 3-08-09 */
    std::fprintf(fp, "\n\nRaw Score Moments (based on IRT fitted distributions)\n\n");
    std::fprintf(fp, " Group Form        mean         sd       skew       kurt\n\n");

    for (i=0;i<4;i++) moments[i] = pinall->IRT_Input->NewForm->mtsng[i];
    std::fprintf(fp, "   1     X    %10.5f %10.5f %10.5f %10.5f\n",
      moments[0], moments[1], moments[2], moments[3]);

    for (i=0;i<4;i++) moments[i] = pinall->IRT_Input->NewForm->mtsog[i];
    std::fprintf(fp, "   2     X    %10.5f %10.5f %10.5f %10.5f\n",
      moments[0], moments[1], moments[2], moments[3]);

    for (i=0;i<4;i++) moments[i] = pinall->IRT_Input->NewForm->mtssg[i];
    std::fprintf(fp, "   s     X    %10.5f %10.5f %10.5f %10.5f\n\n",
      moments[0], moments[1], moments[2], moments[3]);

    for (i=0;i<4;i++) moments[i] = pinall->IRT_Input->OldForm->mtsng[i];
    std::fprintf(fp, "   1     Y    %10.5f %10.5f %10.5f %10.5f\n",
      moments[0], moments[1], moments[2], moments[3]);

    for (i=0;i<4;i++) moments[i] = pinall->IRT_Input->OldForm->mtsog[i];
    std::fprintf(fp, "   2     Y    %10.5f %10.5f %10.5f %10.5f\n",
      moments[0], moments[1], moments[2], moments[3]);

    for (i=0;i<4;i++) moments[i] = pinall->IRT_Input->OldForm->mtssg[i];
    std::fprintf(fp, "   s     Y    %10.5f %10.5f %10.5f %10.5f\n\n",
      moments[0], moments[1], moments[2], moments[3]);

    for(i=1;i<=57;i++) std::fprintf(fp,"-");
  }

  return;

}

/********************************************************************************/

void Print_ESS_QD(FILE *fp, char tt[], struct PDATA *pinall, struct ESS_RESULTS *s,
                  struct RawFitDist *NewForm, struct RawFitDist *OldForm)
/*------------------------------------------------------------------------------*/
{ 
  double moments[4],
         *yctr;                /* conversion table in terms of rounded scale scores */ 
  int i,
      index;         /* location of IRT observed score equivalents in essu and essr */

  /* minimal error checking */

  if(pinall->IRT_Input->method == 'T')
    runerror("This function should not be called when method=='T'"); 

  /* It is assumed here that OldForm->nraws = nscores(maxp,minp,incp), where 
     maxp and minp are defined in comments for ReadSSConvTableForY().
     This assumption should always be true. */

  if(OldForm->nRaws != nscores(pinall->maxp,pinall->minp,pinall->incp)){
    std::printf("OldForm->nRaws = %d but\n" 
             "nscores(pinall->maxp,pinall->minp,pinall->incp) = %d\n" 
             "They should be equal.  One possible explanation is\n"
             "that Wrapper_ESS() was not called previously",
             OldForm->nRaws,nscores(pinall->maxp,pinall->minp,pinall->incp));
    std::printf("\n\nType return to exit\n"); 
    std::getchar();
    std::exit(EXIT_FAILURE);
  }
  
  /* assign value to index */

  if(pinall->IRT_Input->method == 'O') index = 0;   
  else if(pinall->IRT_Input->method == 'A') index = 1;
  else index = -1;                         /* no IRT observed score equating */

  /*** UNROUNDED SCALE SCORE MOMENTS ***/
                                      /* change made in next line on 3-08-09 */
  std::fprintf(fp, "\n\nScale Score Moments (based on IRT fitted distributions)\n"
                  "   TX means based on Form X IRT true score equating\n"
                  "   OX means based on Form X IRT observed score equating\n\n");

  std::fprintf(fp, " Group Form        mean         sd       skew       kurt\n\n");
  std::fprintf(fp, " Unrounded:\n\n"); 

  /* Unrounded scale score moments for old group based on old-form conversion  
     table and old-group quadrature distribution.  These moments might be viewed
     (cautiously) as references for the subsequent moments for IRT true-score and
     observed-score equating. 

     Note that pinall->yct[1][1] = old-form scale score associated with minp.  
  */

  MomentsFromRFD(0,OldForm->nRaws-1,1, &(pinall->yct[1][1]), OldForm->oldFits, moments);
  std::fprintf(fp,"   2     Y    %10.5f %10.5f %10.5f %10.5f\n",
             moments[0], moments[1], moments[2], moments[3]);

  /* Unrounded scale score moments for new-form IRT true score equivalents
     based on new-group quadrature distribution */

  if(pinall->IRT_Input->method == 'T' || pinall->IRT_Input->method == 'A'){
    MomentsFromRFD(0,NewForm->nRaws-1,1, s->essu[0], NewForm->newFits, moments);
    std::fprintf(fp,"   1    TX    %10.5f %10.5f %10.5f %10.5f\n",
               moments[0], moments[1], moments[2], moments[3]);
  }

  /* Unrounded scale score moments for new-form IRT observed score equivalents
     based on new-group quadrature distribution */

  if(index >= 0){
    MomentsFromRFD(0,NewForm->nRaws-1,1, s->essu[index], NewForm->newFits, moments);
    std::fprintf(fp,"   1    OX    %10.5f %10.5f %10.5f %10.5f",
             moments[0], moments[1], moments[2], moments[3]);
  }

  /*** ROUNDED SCALE SCORE MOMENTS ***/

  std::fprintf(fp, "\n\n Rounded:\n\n");

  /* Rounded scale score moments for old group based on old-form conversion  
     table and old-group quadrature distribution.  These moments might be viewed
     (cautiously) as references for the subsequent moments for IRT true-score and
     observed-score equating. */

  yctr = dvector(0,OldForm->nRaws-1);

  for(i=0;i<OldForm->nRaws;i++){    /* round the scale scores in the conv table */
       yctr[i] = std::pow(10.0, static_cast<double>(pinall->round - 1))*
                 (static_cast<int>(pinall->yct[1][i+1]/
                        std::pow(10.0, static_cast<double>(pinall->round - 1)) +.5));
       if(yctr[i] < pinall->lprss) yctr[i] = pinall->lprss;
       if(yctr[i] > pinall->hprss) yctr[i] = pinall->hprss;
  }

  MomentsFromRFD(0,OldForm->nRaws-1,1, yctr, OldForm->oldFits, moments);
  std::fprintf(fp,"   2     Y    %10.5f %10.5f %10.5f %10.5f\n",
             moments[0], moments[1], moments[2], moments[3]);

  /* Rounded scale score moments for new-form IRT true score equivalents
     based on new-group quadrature distribution */

  if(pinall->IRT_Input->method == 'T' || pinall->IRT_Input->method == 'A'){
    MomentsFromRFD(0,NewForm->nRaws-1,1, s->essr[0], NewForm->newFits, moments);
    std::fprintf(fp,"   1    TX    %10.5f %10.5f %10.5f %10.5f\n",
               moments[0], moments[1], moments[2], moments[3]);
  }

  /* Rounded scale score moments for new-form IRT observed score equivalents
     based on new-group quadrature distribution */

  if(index >= 0){
    MomentsFromRFD(0,NewForm->nRaws-1,1, s->essr[index], NewForm->newFits, moments);
    std::fprintf(fp,"   1    OX    %10.5f %10.5f %10.5f %10.5f\n\n",
             moments[0], moments[1], moments[2], moments[3]);
  }

  free_dvector(yctr,0,OldForm->nRaws-1);
  for(i=1;i<=63;i++) std::fprintf(fp,"*");
}
/*******************************************************************************/